#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

int main(void) {
    gmp_randstate_t state;
    mpz_t a, b;
    FILE *fa, *fb;

    // Инициализация чисел
    mpz_init(a);
    mpz_init(b);

    // Инициализация генератора случайных чисел
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));

    // Генерация случайных чисел ~1 КБ (8192 бит)
    mpz_urandomb(a, state, 8192);
    mpz_urandomb(b, state, 8192);

    // Открываем файлы для записи
    fa = fopen("numA.txt", "w");
    fb = fopen("numB.txt", "w");

    if (!fa || !fb) {
        perror("Не удалось открыть файл");
        return 1;
    }

    // Записываем числа в десятичной форме
    mpz_out_str(fa, 10, a); fprintf(fa, "\n");
    mpz_out_str(fb, 10, b); fprintf(fb, "\n");

    fclose(fa);
    fclose(fb);

    printf("Числа сгенерированы и сохранены в numA.txt и numB.txt\n");

    // Очистка
    mpz_clear(a);
    mpz_clear(b);
    gmp_randclear(state);

    return 0;
}
