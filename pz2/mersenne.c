#define _GNU_SOURCE 1
#include "../pz1/big_int.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>

/*
 * Программа для исследования простоты чисел Мерсенна.
 * Реализованы следующие этапы:
 *   1) Построение числа M_p = 2^p - 1 с использованием big_int.c.
 *   2) Набор «стандартных» делимых проверок (чётность, суммы цифр, делимость на 6).
 *   3) Деление на простые числа, полученные с помощью решета Эратосфена.
 *   4) Тест Лука–Лемера для проверки простоты числа Мерсенна.
 * Все комментарии и сообщения на русском языке по требованию задания.
 */

/* Внутренние утилиты ------------------------------------------------------ */

/**
 * Преобразование BigInt в десятичную строку.
 * Возвращает динамически выделенную строку, которую нужно освободить вызовом free().
 */
static char *bi_to_decimal_string(const BigInt *value) {
    char *buffer = NULL;
    size_t size = 0;
    FILE *mem = open_memstream(&buffer, &size);
    if (!mem) {
        return NULL;
    }
    if (!bi_write_to_stream(value, mem, 0)) {
        fclose(mem);
        free(buffer);
        return NULL;
    }
    fclose(mem);
    return buffer;
}

/** Подсчёт суммы десятичных цифр числа. */
static unsigned int decimal_digit_sum(const char *decimal_representation) {
    unsigned int sum = 0;
    for (const char *p = decimal_representation; *p; ++p) {
        if (*p >= '0' && *p <= '9') {
            sum += (unsigned int)(*p - '0');
        }
    }
    return sum;
}

/** Простое чтение беззнакового int из argv. */
static int parse_uint(const char *text, unsigned int *out) {
    char *end = NULL;
    unsigned long value = strtoul(text, &end, 10);
    if (!text || *text == '\0' || (end && *end != '\0')) {
        return 0;
    }
    if (value > UINT_MAX) {
        return 0;
    }
    *out = (unsigned int)value;
    return 1;
}

/** Быстрая проверка простоты для 32-битного показателя. */
static bool is_prime_uint(unsigned int value) {
    if (value < 2) return false;
    if (value == 2) return true;
    if ((value & 1u) == 0u) return false;
    for (unsigned int d = 3; (unsigned long)d * d <= value; d += 2) {
        if (value % d == 0u) {
            return false;
        }
    }
    return true;
}

/** Возвращает наименьший делитель > 1 (для составных показателей). */
static unsigned int smallest_divisor(unsigned int value) {
    if (value % 2u == 0u) return 2u;
    for (unsigned int d = 3; (unsigned long)d * d <= value; d += 2) {
        if (value % d == 0u) {
            return d;
        }
    }
    return value;
}

/**
 * Формирование числа Мерсенна по показателю p.
 * На выходе: res = 2^p - 1. Предполагается, что res заранее инициализирован.
 */
static int build_mersenne(unsigned int p, BigInt *res) {
    BigInt two, exponent, power, one;
    bi_init(&two);
    bi_init(&exponent);
    bi_init(&power);
    bi_init(&one);

    bi_set_u32(&two, 2);
    bi_set_u32(&exponent, p);
    bi_set_u32(&one, 1);

    int ok = bi_pow_bigexp(&two, &exponent, &power);
    if (!ok) {
        bi_free(&one);
        bi_free(&power);
        bi_free(&exponent);
        bi_free(&two);
        return 0;
    }

    bi_sub(&power, &one, res);

    bi_free(&one);
    bi_free(&power);
    bi_free(&exponent);
    bi_free(&two);
    return 1;
}

/* Стандартные проверки ---------------------------------------------------- */

/**
 * Выполнение простых делимых проверок.
 * На экран выводятся поясняющие сообщения, а также вывод о том, прошли ли проверки.
 */
static void run_basic_checks(const BigInt *value, size_t digit_count, unsigned int digit_sum) {
    printf("\nСтандартные проверки на простоту:\n");

    uint32_t mod2 = bi_mod_u32(value, 2);
    printf("  • Проверка на чётность: %s\n", (mod2 == 0) ? "число чётное → составное" : "число нечётное");

    uint32_t mod3 = bi_mod_u32(value, 3);
    uint32_t mod9 = bi_mod_u32(value, 9);
    printf("  • Сумма цифр = %u → %s делится на 3, %s делится на 9\n",
           digit_sum,
           (mod3 == 0) ? "" : "не",
           (mod9 == 0) ? "" : "не");

    uint32_t mod6 = bi_mod_u32(value, 6);
    printf("  • Проверка на 6: остаток %u → %s\n", mod6,
           (mod6 == 1 || mod6 == 5) ? "возможный кандидат (6k±1)" : "число точно составное");

    printf("  • Число имеет %zu десятичных цифр.\n", digit_count);
}

/* Решето Эратосфена ------------------------------------------------------- */

/**
 * Поиск делителей среди первых простых чисел, полученных решетом Эратосфена.
 * Возвращает 1, если делитель не найден до заданной границы, и 0 при обнаружении делителя.
 * В случае нахождения делителя записывает его значение в *divisor.
 */
static int sieve_division_test(const BigInt *value, unsigned int limit, unsigned int *divisor) {
    if (limit < 2) {
        return 1;
    }

    unsigned char *is_prime = (unsigned char*)malloc(limit + 1);
    if (!is_prime) {
        fprintf(stderr, "Не удалось выделить память для решета Эратосфена.\n");
        return 1;
    }
    memset(is_prime, 1, limit + 1);
    is_prime[0] = is_prime[1] = 0;

    for (unsigned int i = 2; (unsigned long long)i * i <= limit; ++i) {
        if (!is_prime[i]) continue;
        for (unsigned int j = i * i; j <= limit; j += i) {
            is_prime[j] = 0;
        }
    }

    int result = 1;
    for (unsigned int p = 2; p <= limit; ++p) {
        if (!is_prime[p]) continue;
        if (bi_mod_u32(value, p) == 0) {
            result = 0;
            if (divisor) *divisor = p;
            break;
        }
    }

    free(is_prime);
    return result;
}

/* Решето Аткина ----------------------------------------------------------- */

/**
 * Поиск делителей с помощью решета Аткина.
 * Возвращает 1, если делитель не найден до заданной границы, и 0 при обнаружении делителя.
 * В случае нахождения делителя записывает его значение в *divisor.
 */
static int sieve_atkin_division_test(const BigInt *value, unsigned int limit, unsigned int *divisor) {
    if (limit < 2) {
        return 1;
    }

    unsigned char *is_prime = (unsigned char*)calloc(limit + 1, sizeof(unsigned char));
    if (!is_prime) {
        fprintf(stderr, "Не удалось выделить память для решета Аткина.\n");
        return 1;
    }

    unsigned int sqrt_limit = (unsigned int)sqrt((double)limit);

    for (unsigned int x = 1; x <= sqrt_limit; ++x) {
        unsigned long long x2 = (unsigned long long)x * (unsigned long long)x;
        for (unsigned int y = 1; y <= sqrt_limit; ++y) {
            unsigned long long y2 = (unsigned long long)y * (unsigned long long)y;

            unsigned long long n = 4ull * x2 + y2;
            if (n <= limit && (n % 12ull == 1ull || n % 12ull == 5ull)) {
                is_prime[n] ^= 1u;
            }

            n = 3ull * x2 + y2;
            if (n <= limit && (n % 12ull == 7ull)) {
                is_prime[n] ^= 1u;
            }

            if (x > y) {
                n = 3ull * x2 - y2;
                if (n <= limit && (n % 12ull == 11ull)) {
                    is_prime[n] ^= 1u;
                }
            }
        }
    }

    if (limit >= 2) is_prime[2] = 1;
    if (limit >= 3) is_prime[3] = 1;

    for (unsigned int r = 5; r <= sqrt_limit; ++r) {
        if (!is_prime[r]) continue;
        unsigned long long r2 = (unsigned long long)r * (unsigned long long)r;
        for (unsigned long long k = r2; k <= limit; k += r2) {
            is_prime[k] = 0;
        }
    }

    int result = 1;
    if (bi_mod_u32(value, 2) == 0) {
        result = 0;
        if (divisor) *divisor = 2;
    } else if (limit >= 3 && bi_mod_u32(value, 3) == 0) {
        result = 0;
        if (divisor) *divisor = 3;
    } else {
        for (unsigned int p = 5; p <= limit; ++p) {
            if (!is_prime[p]) continue;
            if (bi_mod_u32(value, p) == 0) {
                result = 0;
                if (divisor) *divisor = p;
                break;
            }
        }
    }

    free(is_prime);
    return result;
}


/* Тест Лука–Лемера -------------------------------------------------------- */

/**
 * Реализация теста Лука–Лемера для числа Мерсенна M_p.
 * Возвращает 1, если число вероятно простое, и 0 при обнаружении составности.
 */
static int lucas_lehmer_test(unsigned int p, const BigInt *mersenne) {
    if (p < 2) {
        return 0;
    }
    if (p == 2) {
        return 1;
    }

    BigInt s, two, tmp, mod;
    bi_init(&s);
    bi_init(&two);
    bi_init(&tmp);
    bi_init(&mod);

    bi_set_u32(&s, 4);
    bi_set_u32(&two, 2);

    for (unsigned int i = 0; i < p - 2; ++i) {
        BigInt sq;
        bi_init(&sq);
        bi_mul(&s, &s, &sq);       /* s = s^2 */
        bi_sub(&sq, &two, &tmp);   /* tmp = s^2 - 2 */
        bi_mod(&tmp, mersenne, &mod); /* mod = tmp mod M */
        bi_copy(&mod, &s);         /* s = mod */
        bi_free(&sq);
    }

    int is_zero = bi_is_zero(&s);

    bi_free(&mod);
    bi_free(&tmp);
    bi_free(&two);
    bi_free(&s);

    return is_zero;
}

/* Основная программа ------------------------------------------------------ */

int main(int argc, char **argv) {
    unsigned int p = 9941;              /* показатель степени по умолчанию */
    unsigned int sieve_limit = 200000;  /* граница решета Эратосфена */
    const char *output_path = "mersenne.txt"; /* путь для сохранения числа */

    if (argc >= 2) {
        if (!parse_uint(argv[1], &p)) {
            fprintf(stderr, "Ошибка: не удалось разобрать показатель степени '%s'.\n", argv[1]);
            return 1;
        }
    }
    if (argc >= 3) {
        if (!parse_uint(argv[2], &sieve_limit)) {
            fprintf(stderr, "Ошибка: не удалось разобрать границу решета '%s'.\n", argv[2]);
            return 1;
        }
    }
    if (argc >= 4) {
        output_path = argv[3];
        }
    printf("Исследуется число Мерсенна для p = %u.\n", p);

    BigInt mersenne;
    bi_init(&mersenne);

    clock_t start = clock();
    if (!build_mersenne(p, &mersenne)) {
        fprintf(stderr, "Не удалось построить число Мерсенна.\n");
        bi_free(&mersenne);
        return 1;
    }
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;

    char *decimal = bi_to_decimal_string(&mersenne);
    if (!decimal) {
        fprintf(stderr, "Не удалось сформировать десятичное представление числа.\n");
        bi_free(&mersenne);
        return 1;
    }
    size_t digit_count = strlen(decimal);
    unsigned int digit_sum = decimal_digit_sum(decimal);

    printf("Построение числа завершено за %.3f с. Количество цифр: %zu.\n", seconds, digit_count);

    run_basic_checks(&mersenne, digit_count, digit_sum);

    if (!bi_write_to_file(&mersenne, output_path, 1, 0)) {
        fprintf(stderr, "Не удалось записать число в файл '%s'.\n", output_path);
    } else {
        printf("\nЧисло Мерсенна сохранено в файл '%s'.\n", output_path);
    }


    printf("\nРешето Эратосфена до %u:\n", sieve_limit);
    unsigned int divisor = 0;
    start = clock();
    int sieve_result = sieve_division_test(&mersenne, sieve_limit, &divisor);
    end = clock();
    seconds = (double)(end - start) / CLOCKS_PER_SEC;
    if (sieve_result) {
        printf("  Делители до %u не найдены (время %.3f с).\n", sieve_limit, seconds);
        printf("  Для окончательного вывода требуется проверка до √M ≈ 10^{%.0f};\n", p * log10(2.0) / 2.0);
        printf("  при больших p такую границу достичь стандартным перебором крайне трудно.\n");
    } else {
        printf("  Найден делитель %u (время %.3f с) → число составное.\n", divisor, seconds);
    }
        printf("\nРешето Аткина до %u:\n", sieve_limit);
    divisor = 0;
    start = clock();
    sieve_result = sieve_atkin_division_test(&mersenne, sieve_limit, &divisor);
    end = clock();
    seconds = (double)(end - start) / CLOCKS_PER_SEC;
    if (sieve_result) {
        printf("  Делители до %u не найдены (время %.3f с).\n", sieve_limit, seconds);
    } else {
        printf("  Найден делитель %u (время %.3f с) → число составное.\n", divisor, seconds);
    }
    bool prime_exponent = is_prime_uint(p);
    if (!prime_exponent) {
        unsigned int d = smallest_divisor(p);
        printf("\nПоказатель p составной: %u = %u × %u.\n", p, d, p / d);
        printf("По свойствам чисел Мерсенна M_%u гарантированно составно, так как 2^%u − 1 делится на 2^%u − 1.\n",
               p, p, d);
    }
    printf("\nТест Лука–Лемера:\n");
    if (!prime_exponent) {
        printf("  Тест неприменим: показатель p должен быть простым.\n");

    } else {
        start = clock();
        int is_prime = lucas_lehmer_test(p, &mersenne);
        end = clock();
        seconds = (double)(end - start) / CLOCKS_PER_SEC;
        if (is_prime) {
            printf("  Тест Лука–Лемера подтвердил простоту M_%u (%.3f с).\n", p, seconds);
        } else {
            printf("  Тест Лука–Лемера обнаружил составность M_%u (%.3f с).\n", p, seconds);
            printf("  Напоминание: даже при простом p число Мерсенна может оказаться составным (пример: M_11 = 23 × 89).\n");
        }
    }


    free(decimal);
    bi_free(&mersenne);
    return 0;
}
