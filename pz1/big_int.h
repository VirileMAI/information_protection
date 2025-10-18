#ifndef BIG_INT_H
#define BIG_INT_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int sign;          /* 0, +1, -1 */
    uint32_t *d;       /* слова в базе 1e9, младшие сначала */
    size_t n;          /* фактическое кол-во слов */
    size_t cap;        /* ёмкость массива */
} BigInt;

/* Жизненный цикл */
void bi_init(BigInt *a);
void bi_free(BigInt *a);
void bi_set_zero(BigInt *a);

static inline void bi_copy(const BigInt *src, BigInt *dst) {
    if (src == dst) {
        return;
    }

    if (dst->cap < src->n) {
        size_t cap = dst->cap ? dst->cap : 1u;
        while (cap < src->n) {
            cap <<= 1u;
        }
        uint32_t *p = (uint32_t *)realloc(dst->d, cap * sizeof(uint32_t));
        if (!p) {
            fprintf(stderr, "Out of memory in bi_copy\n");
            exit(1);
        }
        dst->d = p;
        dst->cap = cap;
    }

    if (src->n) {
        memcpy(dst->d, src->d, src->n * sizeof(uint32_t));
    }
    dst->n = src->n;
    dst->sign = src->sign;
}

static inline int bi_is_zero(const BigInt *a) {
    return (a->sign == 0) || (a->n == 0);
}

static inline int bi_is_one(const BigInt *a) {
    return (a->sign > 0) && (a->n == 1) && (a->d[0] == 1u);
}

/* Ввод/вывод (десятичный текст) */
int  bi_from_string(BigInt *a, const char *s);
void bi_print(const BigInt *a);
int  bi_read_from_stream(BigInt *a, FILE *fp);
int  bi_write_to_stream(const BigInt *a, FILE *fp, int with_newline);
int  bi_read_from_file(BigInt *a, const char *path);
int  bi_write_to_file(const BigInt *a, const char *path, int with_newline, int append);

/* Сравнение */
int  bi_cmp_abs(const BigInt *a, const BigInt *b); /* -1/0/+1 для |a| ? |b| */
int  bi_cmp(const BigInt *a, const BigInt *b);     /* -1/0/+1 для a ? b */

/* Арифметика */
void bi_add(const BigInt *a, const BigInt *b, BigInt *res);
void bi_sub(const BigInt *a, const BigInt *b, BigInt *res);
void bi_mul(const BigInt *a, const BigInt *b, BigInt *res);
/* Деление с усечением к нулю: res = a / b; возвращает 1 при успехе, 0 если b==0 */
int  bi_div(const BigInt *a, const BigInt *b, BigInt *res);

static inline int bi_mod(const BigInt *a, const BigInt *m, BigInt *res) {
    if (bi_is_zero(m)) {
        return 0;
    }

    BigInt q; bi_init(&q);
    BigInt prod; bi_init(&prod);
    BigInt rem; bi_init(&rem);

    if (!bi_div(a, m, &q)) {
        bi_free(&q);
        bi_free(&prod);
        bi_free(&rem);
        return 0;
    }

    bi_mul(&q, m, &prod);
    bi_sub(a, &prod, &rem);

    if (rem.sign < 0) {
        BigInt mod_abs; bi_init(&mod_abs);
        bi_copy(m, &mod_abs);
        mod_abs.sign = (mod_abs.n ? +1 : 0);

        BigInt tmp; bi_init(&tmp);
        bi_add(&rem, &mod_abs, &tmp);
        bi_copy(&tmp, res);
        bi_free(&tmp);
        bi_free(&mod_abs);
    } else {
        bi_copy(&rem, res);
    }

    if (res->n == 0) {
        res->sign = 0;
    } else if (res->sign < 0) {
        res->sign = +1;
    }

    bi_free(&q);
    bi_free(&prod);
    bi_free(&rem);
    return 1;
}

/* Степень: res = a^e, e>=0 (большой показатель) */
int  bi_pow_bigexp(const BigInt *a, const BigInt *e, BigInt *res);

/* Теория чисел */
void bi_gcd(const BigInt *a, const BigInt *b, BigInt *res);
void bi_lcm(const BigInt *a, const BigInt *b, BigInt *res);

#endif /* BIG_INT_H */
