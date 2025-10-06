#include "big_int.h"
#include <stdio.h>
#include <string.h>

/* читаем два (или три) числа: либо из файлов argv[1], argv[2], argv[3], либо из stdin */
static int read_two(int argc, char **argv, BigInt *a, BigInt *b) {
    if (argc >= 3) {
        if (!bi_read_from_file(a, argv[1])) {
            fprintf(stderr, "Ошибка: не удалось прочитать число из файла '%s'\n", argv[1]);
            return 0;
        }
        if (!bi_read_from_file(b, argv[2])) {
            fprintf(stderr, "Ошибка: не удалось прочитать число из файла '%s'\n", argv[2]);
            return 0;
        }
        return 1;
    } else {
        fprintf(stderr, "Введите два числа (через пробел или перевод строки):\n");
        if (!bi_read_from_stream(a, stdin)) {
            fprintf(stderr, "Ошибка: не удалось прочитать первое число из stdin\n");
            return 0;
        }
        if (!bi_read_from_stream(b, stdin)) {
            fprintf(stderr, "Ошибка: не удалось прочитать второе число из stdin\n");
            return 0;
        }
        return 1;
    }
}

static int read_exp(int argc, char **argv, BigInt *e) {
    if (argc >= 4) {
        if (!bi_read_from_file(e, argv[3])) {
            fprintf(stderr, "Ошибка: не удалось прочитать показатель степени из файла '%s'\n", argv[3]);
            return 0;
        }
        return 1;
    } else {
        printf("Введите показатель степени (необязательно, ENTER чтобы пропустить):\n");
        /* Попытаемся прочитать: если сразу EOF/перевод строки — считаем, что пропущено */
        int c = fgetc(stdin);
        if (c == EOF || c == '\n') {
            return 2; /* код: пропущено пользователем */
        }
        ungetc(c, stdin);
        if (!bi_read_from_stream(e, stdin)) {
            fprintf(stderr, "Ошибка: не удалось прочитать показатель степени\n");
            return 0;
        }
        return 1;
    }
}

int main(int argc, char **argv) {
    BigInt a, b;
    BigInt sum, diff, mul, q, g, l;
    BigInt exp, pwr;
    int have_exp = 0;

    bi_init(&a); bi_init(&b);
    bi_init(&sum); bi_init(&diff); bi_init(&mul); bi_init(&q);
    bi_init(&g); bi_init(&l);
    bi_init(&exp); bi_init(&pwr);

    if (!read_two(argc, argv, &a, &b)) {
        goto cleanup;
    }

    /* пробуем прочитать показатель степени (необязательно) */
    {
        int rc = read_exp(argc, argv, &exp);
        if (rc == 0) goto cleanup;           /* ошибка */
        if (rc == 1) have_exp = 1;           /* есть показатель */
        else have_exp = 0;                   /* пропущено */
    }

    printf("Числа прочитаны успешно.\n");

    /* Сравнение */
    int c = bi_cmp(&a, &b);
    printf("Сравнение: ");
    if (c < 0)      printf("a < b\n");
    else if (c == 0)printf("a == b\n");
    else            printf("a > b\n");

    /* Арифметика */
    bi_add(&a, &b, &sum);
    bi_sub(&a, &b, &diff);
    bi_mul(&a, &b, &mul);

    if (!bi_write_to_file(&sum,  "sum.txt",  1, 0)) fprintf(stderr, "Не удалось записать sum.txt\n");
    if (!bi_write_to_file(&diff, "diff.txt", 1, 0)) fprintf(stderr, "Не удалось записать diff.txt\n");
    if (!bi_write_to_file(&mul,  "mul.txt",  1, 0)) fprintf(stderr, "Не удалось записать mul.txt\n");

    /* Деление (частное с усечением к нулю) */
    if (bi_div(&a, &b, &q)) {
        if (!bi_write_to_file(&q, "div.txt", 1, 0)) fprintf(stderr, "Не удалось записать div.txt\n");
    } else {
        printf("a / b = <деление на 0 невозможно>\n");
    }

    /* НОД/НОК */
    bi_gcd(&a, &b, &g);
    bi_lcm(&a, &b, &l);

    if (!bi_write_to_file(&g, "gcd.txt", 1, 0)) fprintf(stderr, "Не удалось записать gcd.txt\n");
    if (!bi_write_to_file(&l, "lcm.txt", 1, 0)) fprintf(stderr, "Не удалось записать lcm.txt\n");

    /* Возведение в степень (если показатель задан) */
    if (have_exp) {
        if (!bi_pow_bigexp(&a, &exp, &pwr)) {
            fprintf(stderr, "bi_pow_bigexp: показатель отрицательный — операция пропущена\n");
        } else {
            if (!bi_write_to_file(&pwr, "pow.txt", 1, 0)) fprintf(stderr, "Не удалось записать pow.txt\n");
        }
    }

cleanup:
    bi_free(&pwr);
    bi_free(&exp);
    bi_free(&l);
    bi_free(&g);
    bi_free(&q);
    bi_free(&mul);
    bi_free(&diff);
    bi_free(&sum);
    bi_free(&b);
    bi_free(&a);
    return 0;
}
