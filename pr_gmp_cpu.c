#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <signal.h>
#include <setjmp.h>

jmp_buf jump_buffer;

void sigsegv_handler(int sig)
{
    longjmp(jump_buffer, 1);
}

void gcd(mpz_t result, const mpz_t a, const mpz_t b)
{
    mpz_gcd(result, a, b);
}

void mod_mul(mpz_t result, const mpz_t a, const mpz_t b, const mpz_t mod)
{
    mpz_mul(result, a, b);
    mpz_mod(result, result, mod);
}

void pollards_rho_parallel(mpz_t factor, const mpz_t n, int threads)
{
    int found = 0;

#pragma omp parallel num_threads(threads) shared(found, factor)
    {
        mpz_t x, y, c, d;
        gmp_randstate_t state;

        int thread_id = omp_get_thread_num();
        gmp_randinit_mt(state);
        gmp_randseed_ui(state, time(NULL) ^ thread_id ^ clock());

        mpz_inits(x, y, c, d, NULL);

        mpz_t sqrt_n;
        mpz_inits(sqrt_n, NULL);
        mpz_sqrt(sqrt_n, n);

        unsigned long long int range_start = (sqrt_n->_mp_d[0] / threads) * thread_id;
        unsigned long long int range_end = (sqrt_n->_mp_d[0] / threads) * (thread_id + 1);

        for (unsigned long long int i = range_start; i < range_end; ++i)
        {
            if (found)
                break;

            mpz_set_ui(x, i);
            mpz_set(y, x);
            mpz_urandomm(c, state, n);
            mpz_set_ui(d, 1);

            while (mpz_cmp_ui(d, 1) == 0 && !found)
            {
                mod_mul(x, x, x, n);
                mpz_add(x, x, c);
                mpz_mod(x, x, n);

                mod_mul(y, y, y, n);
                mpz_add(y, y, c);
                mpz_mod(y, y, n);
                mod_mul(y, y, y, n);
                mpz_add(y, y, c);
                mpz_mod(y, y, n);

                mpz_sub(d, x, y);
                mpz_abs(d, d);
                gcd(d, d, n);
            }

            if (mpz_cmp(d, n) != 0 && mpz_cmp_ui(d, 1) != 0)
            {
#pragma omp critical
                {
                    if (!found)
                    {
                        mpz_set(factor, d);
                        found = 1;
                    }
                }
            }
        }

        mpz_clears(x, y, c, d, sqrt_n, NULL);
        gmp_randclear(state);
    }
}

int main()
{
    FILE *input_file = fopen("input.txt", "r");
    FILE *output_file = fopen("output_cpu.txt", "w");

    if (!input_file || !output_file)
    {
        printf("Error opening file.\n");
        return 1;
    }

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);

    char n_str[1024];
    int threads = 4;

    signal(SIGSEGV, sigsegv_handler);

    while (fscanf(input_file, "%1023s", n_str) != EOF)
    {
        mpz_set_str(n, n_str, 10);

        fprintf(output_file, "Factoring %s using %d threads...\n", mpz_get_str(NULL, 10, n), threads);

        if (setjmp(jump_buffer) == 0)
        {
            clock_t start = clock();
            pollards_rho_parallel(factor, n, threads);
            clock_t end = clock();

            double time_taken = (double)(end - start) / CLOCKS_PER_SEC;

            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0)
            {
                gmp_fprintf(output_file, "Non-trivial factor: %Zd\n", factor);
            }
            else
            {
                fprintf(output_file, "Failed to find a non-trivial factor.\n");
            }

            fprintf(output_file, "Time taken: %.6f seconds\n\n", time_taken);
        }
        else
        {
            fprintf(output_file, "Segmentation fault encountered, skipping this factorization.\n\n");
        }
    }

    fclose(input_file);
    fclose(output_file);

    mpz_clears(n, factor, NULL);
    return 0;
}
