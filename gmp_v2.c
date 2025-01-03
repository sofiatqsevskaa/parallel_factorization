#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

void gcd(mpz_t result, const mpz_t a, const mpz_t b)
{
    mpz_gcd(result, a, b);
}

void mod_mul(mpz_t result, const mpz_t a, const mpz_t b, const mpz_t mod)
{
    mpz_mul(result, a, b);
    mpz_mod(result, result, mod);
}

void pollards_rho(mpz_t factor, const mpz_t n)
{
    if (mpz_cmp_ui(n, 1) <= 0)
    {
        mpz_set_ui(factor, 1);
        return;
    }
    if (mpz_probab_prime_p(n, 25))
    {
        mpz_set(factor, n);
        return;
    }

    mpz_t x, y, c, d;
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, time(NULL) ^ clock()); // Improved randomness seed

    mpz_inits(x, y, c, d, NULL);

    do
    {
        mpz_urandomm(x, state, n);
        mpz_set(y, x);
        mpz_urandomm(c, state, n);
        mpz_set_ui(d, 1);

        while (mpz_cmp_ui(d, 1) == 0)
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

        if (mpz_cmp(d, n) == 0)
        {
            mpz_set_ui(d, 1);
        }
    } while (mpz_cmp_ui(d, 1) == 0);

    mpz_set(factor, d);

    mpz_clears(x, y, c, d, NULL);
    gmp_randclear(state);
}

int main()
{
    mpz_t n, factor;
    mpz_inits(n, factor, NULL);

    mpz_set_str(n, "151303445596103095416789609238173284269579", 10);

    printf("Factoring %s...\n", mpz_get_str(NULL, 10, n));
    pollards_rho(factor, n);

    if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0)
    {
        gmp_printf("Non-trivial factor: %Zd\n", factor);
    }
    else
    {
        printf("Failed to find a non-trivial factor.\n");
    }

    mpz_clears(n, factor, NULL);
    return 0;
}
