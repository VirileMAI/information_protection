#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

static int file_exists(const char *fname) {
    FILE *f = fopen(fname, "r");
    if (!f) return 0;
    fclose(f);
    return 1;
}

static int read_mpz_from_file(mpz_t z, const char *fname, int base) {
    FILE *f = fopen(fname, "r");
    if (!f) return 0;
    int rc = mpz_inp_str(z, f, base); // rc = число прочитанных цифр; 0 => ошибка
    fclose(f);
    return rc != 0;
}

static int write_mpz_to_file(const mpz_t z, const char *fname, int base) {
    FILE *f = fopen(fname, "w");
    if (!f) return 0;
    if (mpz_out_str(f, base, z) == 0) { fclose(f); return 0; }
    fputc('\n', f);
    fclose(f);
    return 1;
}

// 1 = PASS, 0 = FAIL, -1 = SKIP (нет файла пользователя)
static int check_user_file_eq_mpz(const char *user_fname, const mpz_t ref, int base) {
    if (!file_exists(user_fname)) return -1;
    mpz_t u; mpz_init(u);
    int ok = read_mpz_from_file(u, user_fname, base);
    int pass = ok && (mpz_cmp(u, ref) == 0);
    mpz_clear(u);
    return pass ? 1 : 0;
}

static void print_passfail(const char *label, int status, const char *fname) {
    if (status == -1) {
        printf("[SKIP] %-8s — нет файла %s\n", label, fname);
    } else if (status == 1) {
        printf("[PASS] %-8s — совпадает с %s\n", label, fname);
    } else {
        printf("[FAIL] %-8s — НЕ совпадает с %s\n", label, fname);
    }
}

int main(int argc, char **argv) {
    const char *fa = (argc >= 2) ? argv[1] : "numA.txt";
    const char *fb = (argc >= 3) ? argv[2] : "numB.txt";
    const char *fe = (argc >= 5) ? argv[4] : (file_exists("exp.txt") ? "exp.txt" : NULL); // опционально: 4-й аргумент
    int base = (argc >= 4) ? atoi(argv[3]) : 10;
    if (base < 2 || base > 62) {
        fprintf(stderr, "Основание base=%d недопустимо (2..62)\n", base);
        return 1;
    }

    mpz_t A, B, E, SUM, DIFF, MUL, G, L, Q, R, PWR;
    mpz_inits(A, B, E, SUM, DIFF, MUL, G, L, Q, R, PWR, NULL);

    // 1) Читаем входы
    if (!read_mpz_from_file(A, fa, base)) {
        fprintf(stderr, "Не удалось прочитать число A из '%s' (base=%d)\n", fa, base);
        goto fail;
    }
    if (!read_mpz_from_file(B, fb, base)) {
        fprintf(stderr, "Не удалось прочитать число B из '%s' (base=%d)\n", fb, base);
        goto fail;
    }
    int have_exp = 0;
    if (fe) {
        if (!read_mpz_from_file(E, fe, base)) {
            fprintf(stderr, "Предупреждение: не удалось прочитать показатель степени из '%s' (pow будет пропущен)\n", fe);
        } else {
            have_exp = 1;
        }
    }

    // 2) Эталоны GMP
    mpz_add(SUM, A, B);
    mpz_sub(DIFF, A, B);
    mpz_mul(MUL, A, B);
    mpz_gcd(G, A, B);
    mpz_lcm(L, A, B);

    if (mpz_sgn(B) != 0) {
        mpz_tdiv_qr(Q, R, A, B); // усечение к нулю
    } else {
        mpz_set_ui(Q, 0);
        mpz_set_ui(R, 0);
    }

    // 3) Запись эталонов (кроме pow, он ниже — условный)
    if (!write_mpz_to_file(SUM, "sum.ref.txt",  base) ||
        !write_mpz_to_file(DIFF,"diff.ref.txt", base) ||
        !write_mpz_to_file(MUL, "mul.ref.txt",  base) ||
        !write_mpz_to_file(G,   "gcd.ref.txt",  base) ||
        !write_mpz_to_file(L,   "lcm.ref.txt",  base) ||
        !write_mpz_to_file(Q,   "div_q.ref.txt",base) ||
        !write_mpz_to_file(R,   "div_r.ref.txt",base)) {
        fprintf(stderr, "Ошибка записи эталонных файлов *.ref.txt\n");
        goto fail;
    }

    // 4) Информация по размерам
    size_t bitsA = mpz_sizeinbase(A, 2);
    size_t bitsB = mpz_sizeinbase(B, 2);
    printf("Прочитано:\n");
    printf("  A: %zu бит (~%.2f КБ)\n", bitsA, bitsA/8.0/1024.0);
    printf("  B: %zu бит (~%.2f КБ)\n", bitsB, bitsB/8.0/1024.0);

    // 5) Сравнение стандартных операций
    int st;
    st = check_user_file_eq_mpz("sum.txt",  SUM,  base); print_passfail("sum",  st, "sum.txt");
    st = check_user_file_eq_mpz("diff.txt", DIFF, base); print_passfail("diff", st, "diff.txt");
    st = check_user_file_eq_mpz("mul.txt",  MUL,  base); print_passfail("mul",  st, "mul.txt");
    st = check_user_file_eq_mpz("gcd.txt",  G,    base); print_passfail("gcd",  st, "gcd.txt");
    st = check_user_file_eq_mpz("lcm.txt",  L,    base); print_passfail("lcm",  st, "lcm.txt");

    st = check_user_file_eq_mpz("div.txt",  Q,    base);
    if (st == -1) {
        printf("[SKIP] %-8s — нет файла %s\n", "div(q)", "div.txt");
    } else if (st == 1) {
        printf("[PASS] %-8s — совпадает с %s\n", "div(q)", "div.txt");
    } else {
        printf("[FAIL] %-8s — НЕ совпадает с %s\n", "div(q)", "div.txt");
    }

    // 6) Пауэр: условная проверка
    if (have_exp) {
        if (mpz_sgn(E) < 0) {
            printf("[SKIP] pow     — показатель отрицательный (mpz_pow_ui не поддерживает), пропуск\n");
        } else if (!mpz_fits_ulong_p(E)) {
            printf("[SKIP] pow     — показатель не помещается в unsigned long, пропуск\n");
        } else {
            unsigned long e_ul = mpz_get_ui(E);

            // защитный порог на размер результата (бит): по умолчанию 8 Мбит (~1 МБ)
            const size_t MAX_POW_BITS = 8u * 1024u * 1024u;
            if (e_ul > 0 && (double)bitsA * (double)e_ul > (double)MAX_POW_BITS) {
                printf("[SKIP] pow     — ожидаемый размер результата слишком велик (> ~1 МБ), пропуск\n");
            } else {
                mpz_pow_ui(PWR, A, e_ul);
                if (!write_mpz_to_file(PWR, "pow.ref.txt", base)) {
                    fprintf(stderr, "Ошибка записи pow.ref.txt\n");
                }
                st = check_user_file_eq_mpz("pow.txt", PWR, base);
                print_passfail("pow", st, "pow.txt");
            }
        }
    } else {
        printf("[SKIP] pow     — показатель не задан (нет exp.txt и 4-го аргумента)\n");
    }

    // // 7) Инвариант деления
    // if (mpz_sgn(B) != 0) {
    //     mpz_t T; mpz_init(T);
    //     mpz_mul(T, Q, B);
    //     mpz_add(T, T, R);
    //     int ok = (mpz_cmp(T, A) == 0);
    //     mpz_clear(T);
    //     printf("Инвариант: A = Q*B + R — %s\n", ok ? "OK" : "НАРУШЕН (!!!)");
    // } else {
    //     printf("Деление не выполнялось (B = 0)\n");
    // }

    // printf("\nГотово. Эталоны: *.ref.txt (pow.ref.txt — если расчёт был выполнен)\n");
    mpz_clears(A, B, E, SUM, DIFF, MUL, G, L, Q, R, PWR, NULL);
    return 0;

fail:
    mpz_clears(A, B, E, SUM, DIFF, MUL, G, L, Q, R, PWR, NULL);
    return 1;
}
