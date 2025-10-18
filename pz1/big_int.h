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
void bi_copy(const BigInt *src, BigInt *dst);
int  bi_is_zero(const BigInt *a);
int  bi_is_one(const BigInt *a);

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
int  bi_mod(const BigInt *a, const BigInt *m, BigInt *res);

/* Степень: res = a^e, e>=0 (большой показатель) */
int  bi_pow_bigexp(const BigInt *a, const BigInt *e, BigInt *res);

/* Теория чисел */
void bi_gcd(const BigInt *a, const BigInt *b, BigInt *res);
void bi_lcm(const BigInt *a, const BigInt *b, BigInt *res);

#endif /* BIG_INT_H */
