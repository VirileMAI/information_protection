#include "big_int.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    unsigned passes;
    unsigned fails;
} Stats;

static BigInt CONST_ZERO;
static BigInt CONST_ONE;
static BigInt CONST_TWO;
static int consts_ready = 0;

static void ensure_constants(void) {
    if (consts_ready) {
        return;
    }
    bi_init(&CONST_ZERO);
    bi_init(&CONST_ONE);
    bi_init(&CONST_TWO);
    bi_from_string(&CONST_ZERO, "0");
    bi_from_string(&CONST_ONE, "1");
    bi_from_string(&CONST_TWO, "2");
    consts_ready = 1;
}

static void free_constants(void) {
    if (!consts_ready) {
        return;
    }
    bi_free(&CONST_ZERO);
    bi_free(&CONST_ONE);
    bi_free(&CONST_TWO);
    consts_ready = 0;
}

static void bi_trim_local(BigInt *a) {
    while (a->n > 0 && a->d[a->n - 1] == 0) {
        --a->n;
    }
    if (a->n == 0) {
        a->sign = 0;
    }
}

static void bi_set_ui(BigInt *a, unsigned value) {
    char buffer[16];
    snprintf(buffer, sizeof(buffer), "%u", value);
    bi_from_string(a, buffer);
}

static void bi_set_si(BigInt *a, long long value) {
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%lld", value);
    bi_from_string(a, buffer);
}

static int bi_is_even(const BigInt *a) {
    if (bi_is_zero(a)) {
        return 1;
    }
    return (a->d[0] & 1u) == 0u;
}

static void bi_half(const BigInt *a, BigInt *res) {
    bi_copy(a, res);
    if (bi_is_zero(res)) {
        return;
    }
    unsigned long long carry = 0ULL;
    for (size_t i = res->n; i-- > 0;) {
        unsigned long long cur = res->d[i] + carry * 1000000000ULL;
        res->d[i] = (uint32_t)(cur / 2ULL);
        carry = cur % 2ULL;
    }
    bi_trim_local(res);
    if (!bi_is_zero(res)) {
        res->sign = a->sign;
    }
}

static unsigned bi_mod_small(const BigInt *a, unsigned m) {
    if (m == 0 || bi_is_zero(a)) {
        return 0;
    }
    unsigned long long rem = 0ULL;
    for (size_t i = a->n; i-- > 0;) {
        rem = (rem * 1000000000ULL + a->d[i]) % m;
    }
    if (a->sign < 0 && rem != 0ULL) {
        rem = (unsigned long long)m - rem;
    }
    return (unsigned)rem;
}

static void bi_mod_add(const BigInt *a, const BigInt *b, const BigInt *mod, BigInt *res) {
    BigInt sum;
    bi_init(&sum);
    bi_add(a, b, &sum);
    bi_mod(&sum, mod, res);
    bi_free(&sum);
}

static void bi_mod_sub(const BigInt *a, const BigInt *b, const BigInt *mod, BigInt *res) {
    BigInt diff;
    bi_init(&diff);
    bi_sub(a, b, &diff);
    if (diff.sign < 0) {
        BigInt tmp;
        bi_init(&tmp);
        bi_add(&diff, mod, &tmp);
        bi_copy(&tmp, &diff);
        bi_free(&tmp);
    }
    bi_mod(&diff, mod, res);
    bi_free(&diff);
}

static void bi_mod_mul(const BigInt *a, const BigInt *b, const BigInt *mod, BigInt *res) {
    BigInt prod;
    bi_init(&prod);
    bi_mul(a, b, &prod);
    bi_mod(&prod, mod, res);
    bi_free(&prod);
}

static void bi_mod_mul_small(const BigInt *a, long long k, const BigInt *mod, BigInt *res) {
    BigInt factor;
    bi_init(&factor);
    bi_set_si(&factor, k);
    bi_mod_mul(a, &factor, mod, res);
    bi_free(&factor);
}

static void bi_mod_pow(const BigInt *base, const BigInt *exp, const BigInt *mod, BigInt *res) {
    ensure_constants();
    BigInt result, base_work, exp_copy;
    bi_init(&result);
    bi_init(&base_work);
    bi_init(&exp_copy);
    bi_set_ui(&result, 1);
    bi_mod(base, mod, &base_work);
    bi_copy(exp, &exp_copy);
    while (!bi_is_zero(&exp_copy)) {
        if (!bi_is_even(&exp_copy)) {
            BigInt tmp;
            bi_init(&tmp);
            bi_mod_mul(&result, &base_work, mod, &tmp);
            bi_copy(&tmp, &result);
            bi_free(&tmp);
        }
        BigInt sq;
        bi_init(&sq);
        bi_mod_mul(&base_work, &base_work, mod, &sq);
        bi_copy(&sq, &base_work);
        bi_free(&sq);
        BigInt half;
        bi_init(&half);
        bi_half(&exp_copy, &half);
        bi_copy(&half, &exp_copy);
        bi_free(&half);
    }
    bi_copy(&result, res);
    bi_free(&result);
    bi_free(&base_work);
    bi_free(&exp_copy);
}

static const unsigned MR_BASES[] = {2, 3, 5, 7, 11, 0};

static int miller_rabin_single(const BigInt *n, unsigned base_value) {
    ensure_constants();
    BigInt base;
    bi_init(&base);
    bi_set_ui(&base, base_value);
    if (bi_cmp(n, &CONST_TWO) < 0) {
        bi_free(&base);
        return 0;
    }
    if (bi_cmp(n, &CONST_TWO) == 0) {
        bi_free(&base);
        return 1;
    }
    if (bi_is_even(n)) {
        bi_free(&base);
        return 0;
    }

    BigInt n_minus_one;
    bi_init(&n_minus_one);
    bi_sub(n, &CONST_ONE, &n_minus_one);
    BigInt d;
    bi_init(&d);
    bi_copy(&n_minus_one, &d);
    unsigned s = 0;
    while (bi_is_even(&d)) {
        BigInt half;
        bi_init(&half);
        bi_half(&d, &half);
        bi_copy(&half, &d);
        bi_free(&half);
        ++s;
    }

    BigInt a;
    bi_init(&a);
    bi_mod(&base, n, &a);
    if (bi_cmp(&a, &CONST_TWO) < 0) {
        bi_copy(&CONST_TWO, &a);
    }

    BigInt g;
    bi_init(&g);
    bi_gcd(&a, n, &g);
    if (!bi_is_one(&g)) {
        int equal = (bi_cmp(&g, n) == 0);
        bi_free(&base);
        bi_free(&n_minus_one);
        bi_free(&d);
        bi_free(&a);
        bi_free(&g);
        return equal;
    }

    BigInt x;
    bi_init(&x);
    bi_mod_pow(&a, &d, n, &x);
    if (bi_is_one(&x) || bi_cmp(&x, &n_minus_one) == 0) {
        bi_free(&base);
        bi_free(&n_minus_one);
        bi_free(&d);
        bi_free(&a);
        bi_free(&g);
        bi_free(&x);
        return 1;
    }

    int result = 0;
    for (unsigned r = 1; r < s; ++r) {
        BigInt tmp;
        bi_init(&tmp);
        bi_mod_mul(&x, &x, n, &tmp);
        bi_copy(&tmp, &x);
        bi_free(&tmp);
        if (bi_cmp(&x, &n_minus_one) == 0) {
            result = 1;
            break;
        }
    }

    bi_free(&base);
    bi_free(&n_minus_one);
    bi_free(&d);
    bi_free(&a);
    bi_free(&g);
    bi_free(&x);
    return result;
}

static Stats run_miller_rabin(const BigInt *n, unsigned rounds) {
    Stats stats = {0, 0};
    unsigned idx = 0;
    for (unsigned i = 0; i < rounds; ++i) {
        if (MR_BASES[idx] == 0) {
            idx = 0;
        }
        unsigned base = MR_BASES[idx++];
        if (miller_rabin_single(n, base)) {
            ++stats.passes;
        } else {
            ++stats.fails;
        }
    }
    return stats;
}

static int bi_is_perfect_square(const BigInt *n) {
    ensure_constants();
    if (n->sign < 0) {
        return 0;
    }
    if (bi_is_zero(n) || bi_is_one(n)) {
        return 1;
    }

    BigInt x, y, q, sum;
    bi_init(&x);
    bi_init(&y);
    bi_init(&q);
    bi_init(&sum);

    bi_copy(n, &x);
    BigInt x_plus_one;
    bi_init(&x_plus_one);
    bi_add(&x, &CONST_ONE, &x_plus_one);
    bi_half(&x_plus_one, &y);

    while (bi_cmp(&y, &x) < 0) {
        bi_copy(&y, &x);
        if (!bi_div(n, &x, &q)) {
            bi_free(&x);
            bi_free(&y);
            bi_free(&q);
            bi_free(&sum);
            bi_free(&x_plus_one);
            return 0;
        }
        bi_add(&x, &q, &sum);
        bi_half(&sum, &y);
    }

    BigInt sq;
    bi_init(&sq);
    bi_mul(&x, &x, &sq);
    int is_square = (bi_cmp(&sq, n) == 0);

    bi_free(&x);
    bi_free(&y);
    bi_free(&q);
    bi_free(&sum);
    bi_free(&x_plus_one);
    bi_free(&sq);
    return is_square;
}

static int jacobi_symbol(const BigInt *a_input, const BigInt *n_input) {
    ensure_constants();
    if (bi_is_even(n_input) || n_input->sign <= 0) {
        return 0;
    }

    BigInt a, n;
    bi_init(&a);
    bi_init(&n);
    bi_mod(a_input, n_input, &a);
    bi_copy(n_input, &n);

    int result = 1;
    BigInt half;
    bi_init(&half);

    while (!bi_is_zero(&a)) {
        while (bi_is_even(&a)) {
            bi_half(&a, &half);
            bi_copy(&half, &a);
            unsigned n_mod8 = bi_mod_small(&n, 8);
            if (n_mod8 == 3 || n_mod8 == 5) {
                result = -result;
            }
        }

        if (bi_mod_small(&a, 4) == 3 && bi_mod_small(&n, 4) == 3) {
            result = -result;
        }

        BigInt tmp;
        bi_init(&tmp);
        bi_mod(&n, &a, &tmp);
        bi_copy(&a, &n);
        bi_copy(&tmp, &a);
        bi_free(&tmp);
    }

    int final = (bi_cmp(&n, &CONST_ONE) == 0) ? result : 0;

    bi_free(&a);
    bi_free(&n);
    bi_free(&half);
    return final;
}

static long long select_selfridge_D(const BigInt *n) {
    long long D = 5;
    while (1) {
        BigInt bigD;
        bi_init(&bigD);
        bi_set_si(&bigD, D);

        BigInt absD;
        bi_init(&absD);
        bi_set_si(&absD, (D < 0) ? -D : D);
        BigInt g;
        bi_init(&g);
        bi_gcd(n, &absD, &g);
        if (bi_cmp(&g, &CONST_ONE) > 0 && bi_cmp(&g, n) < 0) {
            bi_free(&bigD);
            bi_free(&absD);
            bi_free(&g);
            return 0;
        }
        bi_free(&absD);
        bi_free(&g);

        int jac = jacobi_symbol(&bigD, n);
        bi_free(&bigD);
        if (jac == -1) {
            return D;
        }
        if (D > 0) {
            D = -D - 2;
        } else {
            D = -D + 2;
        }
    }
}

static int strong_lucas_prp(const BigInt *n) {
    ensure_constants();
    if (bi_cmp(n, &CONST_TWO) < 0) {
        return 0;
    }
    if (bi_cmp(n, &CONST_TWO) == 0) {
        return 1;
    }
    if (bi_is_even(n) || bi_is_perfect_square(n)) {
        return 0;
    }

    long long D = select_selfridge_D(n);
    if (D == 0) {
        return 0;
    }
    long long Q_val = (1 - D) / 4;

    BigInt Q, D_big;
    bi_init(&Q);
    bi_init(&D_big);
    bi_set_si(&Q, Q_val);
    bi_set_si(&D_big, D);

    BigInt n_plus_one, d;
    bi_init(&n_plus_one);
    bi_init(&d);
    bi_add(n, &CONST_ONE, &n_plus_one);
    bi_copy(&n_plus_one, &d);

    unsigned s = 0;
    BigInt half;
    bi_init(&half);
    while (bi_is_even(&d)) {
        bi_half(&d, &half);
        bi_copy(&half, &d);
        ++s;
    }

    BigInt inv2;
    bi_init(&inv2);
    BigInt tmp;
    bi_init(&tmp);
    bi_add(n, &CONST_ONE, &tmp);
    bi_half(&tmp, &inv2);

    // Prepare binary representation of d
    BigInt U, V, Qk;
    bi_init(&U);
    bi_init(&V);
    bi_init(&Qk);
    bi_set_zero(&U);
    bi_set_ui(&V, 2);
    bi_set_ui(&Qk, 1);

    unsigned char *bits = NULL;
    size_t bits_size = 0;
    size_t bits_cap = 0;
    BigInt d_copy;
    bi_init(&d_copy);
    bi_copy(&d, &d_copy);
    while (!bi_is_zero(&d_copy)) {
        if (bits_size == bits_cap) {
            bits_cap = bits_cap ? bits_cap * 2 : 32;
            unsigned char *new_bits = (unsigned char *)realloc(bits, bits_cap);
            if (!new_bits) {
                fprintf(stderr, "Out of memory in strong_lucas_prp\n");
                free(bits);
                bits = NULL;
                bits_cap = 0;
                bi_free(&Q);
                bi_free(&D_big);
                bi_free(&n_plus_one);
                bi_free(&d);
                bi_free(&half);
                bi_free(&inv2);
                bi_free(&tmp);
                bi_free(&d_copy);
                bi_free(&U);
                bi_free(&V);
                bi_free(&Qk);
                return 0;
            }
            bits = new_bits;
        }
        bits[bits_size++] = bi_is_even(&d_copy) ? 0 : 1;
        bi_half(&d_copy, &half);
        bi_copy(&half, &d_copy);
    }

    for (size_t i = bits_size; i-- > 0;) {
        // Square step
        BigInt U_tmp;
        bi_init(&U_tmp);
        bi_mod_mul(&U, &V, n, &U_tmp);

        BigInt V_sq;
        bi_init(&V_sq);
        bi_mod_mul(&V, &V, n, &V_sq);
        BigInt twoQk;
        bi_init(&twoQk);
        bi_mod_mul_small(&Qk, 2, n, &twoQk);
        BigInt V_tmp;
        bi_init(&V_tmp);
        bi_mod_sub(&V_sq, &twoQk, n, &V_tmp);

        BigInt Qk_sq;
        bi_init(&Qk_sq);
        bi_mod_mul(&Qk, &Qk, n, &Qk_sq);

        bi_copy(&U_tmp, &U);
        bi_copy(&V_tmp, &V);
        bi_copy(&Qk_sq, &Qk);

        if (bits[i]) {
            BigInt U_plus_V;
            bi_init(&U_plus_V);
            bi_mod_add(&U, &V, n, &U_plus_V);
            BigInt U_new;
            bi_init(&U_new);
            bi_mod_mul(&U_plus_V, &inv2, n, &U_new);

            BigInt D_mul_U;
            bi_init(&D_mul_U);
            bi_mod_mul(&D_big, &U, n, &D_mul_U);
            BigInt D_mul_U_plus_V;
            bi_init(&D_mul_U_plus_V);
            bi_mod_add(&D_mul_U, &V, n, &D_mul_U_plus_V);
            BigInt V_new;
            bi_init(&V_new);
            bi_mod_mul(&D_mul_U_plus_V, &inv2, n, &V_new);

            BigInt Qk_new;
            bi_init(&Qk_new);
            bi_mod_mul(&Qk, &Q, n, &Qk_new);

            bi_copy(&U_new, &U);
            bi_copy(&V_new, &V);
            bi_copy(&Qk_new, &Qk);

            bi_free(&U_plus_V);
            bi_free(&U_new);
            bi_free(&D_mul_U);
            bi_free(&D_mul_U_plus_V);
            bi_free(&V_new);
            bi_free(&Qk_new);
        }

        bi_free(&U_tmp);
        bi_free(&V_sq);
        bi_free(&twoQk);
        bi_free(&V_tmp);
        bi_free(&Qk_sq);
    }

    int is_probably_prime = 0;
    if (bi_is_zero(&U) || bi_is_zero(&V)) {
        is_probably_prime = 1;
    }

    for (unsigned r = 1; !is_probably_prime && r < s; ++r) {
        BigInt V_sq;
        bi_init(&V_sq);
        bi_mod_mul(&V, &V, n, &V_sq);
        BigInt twoQk;
        bi_init(&twoQk);
        bi_mod_mul_small(&Qk, 2, n, &twoQk);
        BigInt V_next;
        bi_init(&V_next);
        bi_mod_sub(&V_sq, &twoQk, n, &V_next);

        BigInt Qk_sq;
        bi_init(&Qk_sq);
        bi_mod_mul(&Qk, &Qk, n, &Qk_sq);

        bi_copy(&V_next, &V);
        bi_copy(&Qk_sq, &Qk);

        if (bi_is_zero(&V)) {
            is_probably_prime = 1;
        }

        bi_free(&V_sq);
        bi_free(&twoQk);
        bi_free(&V_next);
        bi_free(&Qk_sq);
    }

    bi_free(&Q);
    bi_free(&D_big);
    bi_free(&n_plus_one);
    bi_free(&d);
    bi_free(&half);
    bi_free(&inv2);
    bi_free(&tmp);
    bi_free(&d_copy);
    bi_free(&U);
    bi_free(&V);
    bi_free(&Qk);
    free(bits);

    return is_probably_prime;
}

static int baillie_psw(const BigInt *n) {
    ensure_constants();
    if (bi_cmp(n, &CONST_TWO) < 0) {
        return 0;
    }
    if (bi_cmp(n, &CONST_TWO) == 0) {
        return 1;
    }
    if (bi_is_even(n) || bi_is_perfect_square(n)) {
        return 0;
    }

    if (!miller_rabin_single(n, 2)) {
        return 0;
    }
    return strong_lucas_prp(n);
}

static Stats run_strong_lucas(const BigInt *n, unsigned rounds) {
    Stats stats = {0, 0};
    for (unsigned i = 0; i < rounds; ++i) {
        if (strong_lucas_prp(n)) {
            ++stats.passes;
        } else {
            ++stats.fails;
        }
    }
    return stats;
}

static Stats run_baillie_psw(const BigInt *n, unsigned rounds) {
    Stats stats = {0, 0};
    for (unsigned i = 0; i < rounds; ++i) {
        if (baillie_psw(n)) {
            ++stats.passes;
        } else {
            ++stats.fails;
        }
    }
    return stats;
}

static void print_stats(const char *label, const BigInt *n, const Stats *stats) {
    (void)n;
    printf("  %-15s : passes=%3u fails=%3u\n", label, stats->passes, stats->fails);
}

int main(int argc, char **argv) {
    ensure_constants();
    BigInt prime;
    BigInt composite;
    bi_init(&prime);
    bi_init(&composite);
    bi_from_string(&prime, "32416190071");
    bi_from_string(&composite, "32416190069");

    if (argc > 1) {
        if (!bi_from_string(&prime, argv[1])) {
            fprintf(stderr, "Не удалось прочитать первое число\n");
            return EXIT_FAILURE;
        }
    }
    if (argc > 2) {
        if (!bi_from_string(&composite, argv[2])) {
            fprintf(stderr, "Не удалось прочитать второе число\n");
            return EXIT_FAILURE;
        }
    }

    printf("Проверка вероятно простого: ");
    bi_print(&prime);
    printf("\n");
    Stats prime_mr = run_miller_rabin(&prime, 100);
    Stats prime_lucas = run_strong_lucas(&prime, 100);
    Stats prime_bpsw = run_baillie_psw(&prime, 100);
    print_stats("Miller-Rabin", &prime, &prime_mr);
    print_stats("Lucas", &prime, &prime_lucas);
    print_stats("Baillie-PSW", &prime, &prime_bpsw);
    printf("  Вердикт          : %s\n", (prime_bpsw.fails == 0) ? "вероятно простое" : "составное");

    printf("\nПроверка составного: ");
    bi_print(&composite);
    printf("\n");
    Stats comp_mr = run_miller_rabin(&composite, 100);
    Stats comp_lucas = run_strong_lucas(&composite, 100);
    Stats comp_bpsw = run_baillie_psw(&composite, 100);
    print_stats("Miller-Rabin", &composite, &comp_mr);
    print_stats("Lucas", &composite, &comp_lucas);
    print_stats("Baillie-PSW", &composite, &comp_bpsw);
    printf("  Вердикт          : %s\n", (comp_bpsw.fails == 100) ? "составное" : "возможно простое");

    bi_free(&prime);
    bi_free(&composite);
    free_constants();
    return 0;
}

