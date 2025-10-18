#include "big_int.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#define BI_BASE  1000000000u   /* 10^9 */
#define BI_CHUNK 9

// Вспомогательные функции
static void bi_reserve(BigInt *a, size_t need){
    if (a->cap >= need) return;
    size_t cap = a->cap ? a->cap : 4;
    while (cap < need) cap <<= 1;
    uint32_t *p = (uint32_t*)realloc(a->d, cap*sizeof(uint32_t));
    if(!p){ fprintf(stderr,"OOM\n"); exit(1); }
    a->d = p; a->cap = cap;
}
static void bi_trim(BigInt *a){
    while (a->n && a->d[a->n-1]==0) a->n--;
    if (!a->n) a->sign = 0;
}
static const char* skip_ws(const char* s){
    while(*s && isspace((unsigned char)*s)) ++s;
    return s;
}

/* helper для парсера: res = res*m + add, где m<=1e9, add<1e9 */
static void mul_small_add_inplace(BigInt *a, uint32_t m, uint32_t add){
    if (a->sign==0 && a->n==0) {
        bi_reserve(a,1);
        a->d[0]=0; a->n=1; a->sign=+1;
    }
    uint64_t carry=0;
    for(size_t i=0;i<a->n;++i){
        uint64_t v = (uint64_t)a->d[i]*m + carry;
        a->d[i]=(uint32_t)(v%BI_BASE);
        carry = v/BI_BASE;
    }
    if (carry){ bi_reserve(a, a->n+1); a->d[a->n++]=(uint32_t)carry; }
    /* + add (младшее слово) */
    uint64_t s = (uint64_t)a->d[0] + add;
    a->d[0]=(uint32_t)(s%BI_BASE);
    uint64_t c = s/BI_BASE;
    for(size_t i=1; c && i<a->n; ++i){
        uint64_t z = (uint64_t)a->d[i] + c;
        a->d[i]=(uint32_t)(z%BI_BASE);
        c = z/BI_BASE;
    }
    if (c){ bi_reserve(a, a->n+1); a->d[a->n++]=(uint32_t)c; }
}

/* умножение на слово: out = a * m; out != a допустимо */
static void mul_word(const BigInt *a, uint32_t m, BigInt *out) {
    if (a->sign == 0 || m == 0) { bi_set_zero(out); return; }
    bi_reserve(out, a->n + 1);
    uint64_t carry = 0;
    for (size_t i = 0; i < a->n; ++i) {
        uint64_t cur = (uint64_t)a->d[i] * m + carry;
        out->d[i] = (uint32_t)(cur % BI_BASE);
        carry = cur / BI_BASE;
    }
    out->n = a->n;
    if (carry) out->d[out->n++] = (uint32_t)carry;
    out->sign = a->sign;
    bi_trim(out);
}

//Возведение в квадрат
static void sqr_school(const BigInt *a, BigInt *out) {
    if (a->sign == 0) { bi_set_zero(out); return; }
    size_t n = a->n;
    bi_reserve(out, 2*n + 1);
    size_t out_len = 2*n + 1;
    memset(out->d, 0, out_len * sizeof(uint32_t));

    for (size_t i = 0; i < n; ++i) {
        uint64_t carry = out->d[2*i] + (uint64_t)a->d[i] * (uint64_t)a->d[i];
        out->d[2*i] = (uint32_t)(carry % BI_BASE);
        uint64_t hi = carry / BI_BASE;

        uint64_t c = hi;
        for (size_t j = i + 1; j < n; ++j) {
            __uint128_t cur = (__uint128_t)out->d[i + j]
                            + 2u * (__uint128_t)a->d[i] * (__uint128_t)a->d[j]
                            + c;
            out->d[i + j] = (uint32_t)((uint64_t)cur % BI_BASE);
            c = (uint64_t)cur / BI_BASE;
        }
        size_t pos = i + n;
        while (c) {
            __uint128_t cur = (__uint128_t)out->d[pos] + c;
            out->d[pos] = (uint32_t)((uint64_t)cur % BI_BASE);
            c = (uint64_t)cur / BI_BASE;
            ++pos;
        }
    }

    out->n = out_len;
    out->sign = +1;
    bi_trim(out);
}

//Вычитание по модулю, |a|>=|b|
static void sub_abs_ge(const BigInt *a, const BigInt *b, BigInt *c){
    bi_reserve(c, a->n);
    int64_t carry=0;
    size_t i=0;
    for(; i<a->n; ++i){
        int64_t ai=a->d[i];
        int64_t bi=(i<b->n)?b->d[i]:0;
        int64_t s = ai - bi - carry;
        if (s<0){ s+=BI_BASE; carry=1; } else carry=0;
        c->d[i]=(uint32_t)s;
    }
    c->n=a->n;
    bi_trim(c);
}

//Основа
void bi_init(BigInt *a){ a->sign=0; a->d=NULL; a->n=0; a->cap=0; }
void bi_free(BigInt *a){ free(a->d); a->d=NULL; a->n=a->cap=0; a->sign=0; }
void bi_set_zero(BigInt *a){ a->sign=0; a->n=0; }

//Парсинг и вывод
int bi_from_string(BigInt *a, const char *s){
    if(!s) return 0;
    s = skip_ws(s);
    int neg = (*s=='-');
    if (*s=='-' || *s=='+') ++s;

    const char* p = s;
    int has_digit = 0;
    while (isdigit((unsigned char)*p)){ has_digit=1; ++p; }
    const char* t = skip_ws(p);
    if (*t!='\0') return 0;

    bi_set_zero(a);
    if(!has_digit){ return 1; } /* пусто -> ноль */

    size_t len = (size_t)(p - s);
    size_t full = len / BI_CHUNK;
    size_t head = len % BI_CHUNK;

    /* сначала head (1..8 цифр) */
    if(head){
        uint32_t m=1; for(size_t i=1;i<head;i++) m*=10u;
        uint32_t val=0;
        for(size_t i=0;i<head;i++) { val = val*10u + (uint32_t)(s[i]-'0'); }
        bi_reserve(a,1); a->n=1; a->d[0]=0; a->sign=+1;
        mul_small_add_inplace(a, m, val);
    } else {
        bi_reserve(a,1); a->n=1; a->d[0]=0; a->sign=+1;
    }

    /* затем блоки по 9 цифр: res = res*1e9 + val */
    for(size_t k=0;k<full;k++){
        const char* q = s + head + k*BI_CHUNK;
        uint32_t val=0;
        for(size_t i=0;i<BI_CHUNK;i++) val = val*10u + (uint32_t)(q[i]-'0');
        mul_small_add_inplace(a, 1000000000u, val);
    }

    bi_trim(a);
    if (a->n==0) a->sign=0; else if (neg) a->sign=-1;
    return 1;
}

void bi_print(const BigInt *a){
    if (a->sign==0){ printf("0"); return; }
    if (a->sign<0) putchar('-');
    printf("%u", a->d[a->n-1]);
    for (size_t i=a->n-1; i-- > 0; ){
        printf("%09u", a->d[i]);
    }
}

/* потоковый ввод: читаем в строковый буфер динамически (до EOF/пробела после числа) */
int bi_read_from_stream(BigInt *a, FILE *fp){
    if(!fp) return 0;
    int c;
    do { c=fgetc(fp); if (c==EOF) return 0; } while (isspace(c));
    size_t cap=1024, len=0; char *buf=(char*)malloc(cap);
    if(!buf) return 0;
    if (c=='+'||c=='-'||isdigit(c)) { buf[len++]=(char)c; }
    else { free(buf); return 0; }
    while(1){
        int d=fgetc(fp);
        if (d==EOF || isspace(d)){
            buf[len]='\0';
            int ok=bi_from_string(a, buf);
            free(buf);
            return ok?1:0;
        }
        if (!isdigit(d)){ free(buf); return 0; }
        if (len+1>=cap){ cap*=2; char* nb=(char*)realloc(buf,cap); if(!nb){ free(buf); return 0; } buf=nb; }
        buf[len++]=(char)d;
    }
}

int bi_write_to_stream(const BigInt *a, FILE *fp, int with_newline){
    if(!fp) return 0;
    if (a->sign==0){ fputc('0', fp); if(with_newline) fputc('\n',fp); return 1; }
    if (a->sign<0) fputc('-', fp);
    fprintf(fp, "%u", a->d[a->n-1]);
    for (size_t i=a->n-1; i-- > 0; ){
        fprintf(fp, "%09u", a->d[i]);
    }
    if (with_newline) fputc('\n', fp);
    return 1;
}

int bi_read_from_file(BigInt *a, const char *path){
    FILE *f=fopen(path,"rb"); if(!f) return 0;
    int ok = bi_read_from_stream(a,f);
    fclose(f); return ok;
}
int bi_write_to_file(const BigInt *a, const char *path, int with_newline, int append){
    FILE *f=fopen(path, append?"ab":"wb"); if(!f) return 0;
    int ok = bi_write_to_stream(a,f,with_newline);
    fclose(f); return ok;
}

//Сравнение
int bi_cmp_abs(const BigInt *a, const BigInt *b){
    if (a->n != b->n) return (a->n < b->n) ? -1 : +1;
    for (size_t i=a->n; i-- > 0; ){
        if (a->d[i]!=b->d[i]) return (a->d[i]<b->d[i])?-1:+1;
    }
    return 0;
}
int bi_cmp(const BigInt *a, const BigInt *b){
    if (a->sign != b->sign) return (a->sign < b->sign)?-1:+1;
    if (a->sign==0) return 0;
    int c = bi_cmp_abs(a,b);
    return (a->sign>0)?c:-c;
}

//Вычитание, сложение
void bi_add(const BigInt *a, const BigInt *b, BigInt *c){
    if (a->sign==0){
        bi_reserve(c,b->n); memcpy(c->d,b->d,b->n*sizeof(uint32_t)); c->n=b->n; c->sign=b->sign; return;
    }
    if (b->sign==0){ bi_reserve(c,a->n); memcpy(c->d,a->d,a->n*sizeof(uint32_t)); c->n=a->n; c->sign=a->sign; return; }
    if (a->sign == b->sign){
        size_t n = (a->n>b->n)?a->n:b->n;
        bi_reserve(c, n+1);
        uint64_t carry=0;
        size_t i=0;
        for(; i<n; ++i){
            uint64_t s = carry;
            if(i<a->n) s+=a->d[i];
            if(i<b->n) s+=b->d[i];
            c->d[i]=(uint32_t)(s%BI_BASE);
            carry = s/BI_BASE;
        }
        if(carry) c->d[i++]=(uint32_t)carry;
        c->n=i; c->sign=a->sign; bi_trim(c);
    } else {
        int cmp = bi_cmp_abs(a,b);
        if (cmp==0){ bi_set_zero(c); return; }
        if (cmp>0){ sub_abs_ge(a,b,c); c->sign=a->sign; }
        else { sub_abs_ge(b,a,c); c->sign=b->sign; }
    }
}
void bi_sub(const BigInt *a, const BigInt *b, BigInt *c){
    BigInt nb = *b; nb.sign = -b->sign;
    bi_add(a,&nb,c);
}

//Умножение
void bi_mul(const BigInt *a, const BigInt *b, BigInt *res) {
    if (a->sign == 0 || b->sign == 0) { bi_set_zero(res); return; }

    const BigInt *u = a, *v = b;
    if (a->n < b->n) { u = b; v = a; }

    BigInt tmp, *out = res;
    if (res == a || res == b) { bi_init(&tmp); out = &tmp; }

    if (v->n == 1) {
        mul_word(u, v->d[0], out);
        out->sign = u->sign * v->sign;
        bi_trim(out);
        if (out == &tmp) {
            bi_reserve(res, out->n);
            memcpy(res->d, out->d, out->n * sizeof(uint32_t));
            res->n = out->n; res->sign = out->sign;
            bi_free(&tmp);
        }
        return;
    }

    if (u == v) {
        sqr_school(u, out);
        out->sign = +1;
        if (out == &tmp) {
            bi_reserve(res, out->n);
            memcpy(res->d, out->d, out->n * sizeof(uint32_t));
            res->n = out->n; res->sign = out->sign;
            bi_free(&tmp);
        }
        return;
    }

    size_t n_u = u->n, n_v = v->n;
    bi_reserve(out, n_u + n_v + 1);
    size_t out_len = n_u + n_v + 1;
    memset(out->d, 0, out_len * sizeof(uint32_t));

    for (size_t i = 0; i < n_u; ++i) {
        uint64_t carry = 0;
        for (size_t j = 0; j < n_v; ++j) {
            uint64_t cur = (uint64_t)out->d[i + j]
                         + (uint64_t)u->d[i] * (uint64_t)v->d[j]
                         + carry;
            out->d[i + j] = (uint32_t)(cur % BI_BASE);
            carry         = cur / BI_BASE;
        }
        uint64_t cur2 = (uint64_t)out->d[i + n_v] + carry;
        out->d[i + n_v] = (uint32_t)(cur2 % BI_BASE);
        uint64_t c2 = cur2 / BI_BASE;
        if (c2) {
            out->d[i + n_v + 1] += (uint32_t)c2;
        }
    }

    out->n = out_len;
    out->sign = u->sign * v->sign;
    bi_trim(out);

    if (out == &tmp) {
        bi_reserve(res, out->n);
        memcpy(res->d, out->d, out->n * sizeof(uint32_t));
        res->n = out->n; res->sign = out->sign;
        bi_free(&tmp);
    }
}

//Деление
static void sub_shifted(BigInt* a, const BigInt* b, size_t shift){
    int64_t carry=0;
    for(size_t i=0;i<a->n;++i){
        int64_t bi = (i>=shift && (i-shift)<b->n) ? b->d[i-shift] : 0;
        int64_t s = (int64_t)a->d[i] - bi - carry;
        if (s<0){ s+=BI_BASE; carry=1; } else carry=0;
        a->d[i]=(uint32_t)s;
    }
    bi_trim(a);
}
int bi_div(const BigInt *A, const BigInt *B, BigInt *Q){
    if (bi_is_zero(B)) return 0;
    if (bi_is_zero(A)){ bi_set_zero(Q); return 1; }

    BigInt a; bi_init(&a); bi_copy(A,&a); a.sign=+1;
    BigInt b; bi_init(&b); bi_copy(B,&b); b.sign=+1;

    if (bi_cmp_abs(&a,&b) < 0){ bi_set_zero(Q); bi_free(&a); bi_free(&b); return 1; }

    BigInt q; bi_init(&q);
    size_t shift = (a.n > b.n) ? (a.n - b.n) : 0;
    bi_reserve(&q, shift+1);
    q.n = shift+1; memset(q.d, 0, q.n*sizeof(uint32_t)); q.sign=+1;

    for(size_t k=shift+1; k-- > 0; ){
        uint32_t lo=0, hi=BI_BASE-1, best=0;
        while (lo<=hi){
            uint32_t mid = (lo+hi)/2;

            BigInt t; bi_init(&t);
            if (mid){
                BigInt bm; bi_init(&bm);
                bi_reserve(&bm, b.n+1); memset(bm.d,0,(b.n+1)*sizeof(uint32_t)); bm.n=b.n; bm.sign=+1;
                uint64_t carry=0;
                for(size_t i=0;i<b.n;++i){
                    uint64_t cur = (uint64_t)b.d[i]*mid + carry;
                    bm.d[i]=(uint32_t)(cur%BI_BASE);
                    carry=cur/BI_BASE;
                }
                if (carry) bm.d[bm.n++]=(uint32_t)carry;

                bi_reserve(&t, bm.n + k);
                memset(t.d, 0, (bm.n+k)*sizeof(uint32_t));
                memcpy(t.d + k, bm.d, bm.n*sizeof(uint32_t));
                t.n = bm.n + k; t.sign=+1;
                bi_free(&bm);
            } else {
                bi_reserve(&t,1); t.n=1; t.d[0]=0; t.sign=0;
            }

            int cmp = bi_cmp_abs(&a, &t);
            if (cmp>=0){
                best=mid; lo=mid+1;
                bi_free(&t);
            } else {
                bi_free(&t);
                if (mid==0) break;
                hi=mid-1;
            }
        }

        if (best){
            BigInt bm; bi_init(&bm);
            bi_reserve(&bm, b.n+1); memset(bm.d,0,(b.n+1)*sizeof(uint32_t)); bm.n=b.n; bm.sign=+1;
            uint64_t carry=0;
            for(size_t i=0;i<b.n;++i){
                uint64_t cur = (uint64_t)b.d[i]*best + carry;
                bm.d[i]=(uint32_t)(cur%BI_BASE);
                carry=cur/BI_BASE;
            }
            if (carry) bm.d[bm.n++]=(uint32_t)carry;

            sub_shifted(&a, &bm, k);
            bi_free(&bm);
        }
        q.d[k]=best;
        bi_trim(&a);
    }

    bi_trim(&q);
    q.sign = (A->sign==B->sign)? +1 : -1;
    if (q.n==0) q.sign=0;
    bi_copy(&q, Q);
    bi_free(&q); bi_free(&a); bi_free(&b);
    return 1;
}

//Степень
int bi_pow_bigexp(const BigInt *a, const BigInt *e, BigInt *res){
    if (e->sign<0) return 0;
    bi_reserve(res,1); res->n=1; res->d[0]=1; res->sign=+1;

    BigInt base; bi_init(&base); bi_copy(a,&base); base.sign = (base.sign>=0)? +1 : -1;
    BigInt exp;  bi_init(&exp);  bi_copy(e,&exp);  exp.sign=+1;

    while (exp.n){
        if (exp.d[0] & 1u){
            BigInt tmp; bi_init(&tmp);
            bi_mul(res,&base,&tmp);
            bi_copy(&tmp,res); bi_free(&tmp);
        }
        BigInt sq; bi_init(&sq);
        bi_mul(&base,&base,&sq);
        bi_copy(&sq,&base); bi_free(&sq);

        uint32_t carry=0;
        for (size_t i=exp.n; i-- > 0; ){
            uint64_t cur = ((uint64_t)carry*BI_BASE + exp.d[i]);
            exp.d[i] = (uint32_t)(cur/2);
            carry = (uint32_t)(cur%2);
        }
        bi_trim(&exp);
    }

    if (a->sign<0){
        if (e->n && (e->d[0] & 1u)) res->sign = -res->sign;
    }
    bi_free(&base); bi_free(&exp);
    return 1;
}

//НОД, НОК
static int is_even(const BigInt* x){ return (x->n==0) || ((x->d[0] & 1u)==0); }
static void rshift1(BigInt* x){
    uint32_t carry=0;
    for (size_t i=x->n; i-- > 0; ){
        uint64_t cur = (uint64_t)carry*BI_BASE + x->d[i];
        x->d[i] = (uint32_t)(cur/2);
        carry = (uint32_t)(cur%2);
    }
    bi_trim(x);
}
void bi_gcd(const BigInt *A, const BigInt *B, BigInt *res){
    if (bi_is_zero(A)){ bi_copy(B,res); res->sign = (res->n?+1:0); return; }
    if (bi_is_zero(B)){ bi_copy(A,res); res->sign = (res->n?+1:0); return; }
    BigInt a; bi_init(&a); bi_copy(A,&a); a.sign=+1;
    BigInt b; bi_init(&b); bi_copy(B,&b); b.sign=+1;

    size_t k=0;
    while (is_even(&a) && is_even(&b)){ rshift1(&a); rshift1(&b); ++k; }
    while (is_even(&a)) rshift1(&a);

    while (!bi_is_zero(&b)){
        while (is_even(&b)) rshift1(&b);
        int cmp = bi_cmp_abs(&a,&b);
        if (cmp>0){ BigInt t=a; a=b; b=t; }
        sub_abs_ge(&b,&a,&b);
    }

    if (k==0){ bi_copy(&a,res); res->sign=+1; }
    else {
        bi_copy(&a,res); res->sign=+1;
        for(size_t i=0;i<k;i++){
            uint64_t carry=0;
            for(size_t j=0;j<res->n;++j){
                uint64_t cur = (uint64_t)res->d[j]*2 + carry;
                res->d[j]=(uint32_t)(cur%BI_BASE);
                carry = cur/BI_BASE;
            }
            if (carry){ bi_reserve(res,res->n+1); res->d[res->n++]=(uint32_t)carry; }
        }
    }
    bi_free(&a); bi_free(&b);
}

void bi_lcm(const BigInt *a, const BigInt *b, BigInt *res){
    if (bi_is_zero(a) || bi_is_zero(b)){ bi_set_zero(res); return; }
    BigInt g; bi_init(&g); bi_gcd(a,b,&g);
    BigInt absa; bi_init(&absa); bi_copy(a,&absa); absa.sign=+1;
    BigInt t; bi_init(&t);
    bi_div(&absa, &g, &t);
    BigInt absb; bi_init(&absb); bi_copy(b,&absb); absb.sign=+1;
    bi_mul(&t,&absb,res);
    res->sign=+1;
    bi_free(&absb); bi_free(&t); bi_free(&absa); bi_free(&g);
}
