#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <cuda.h>
#include <time.h>

__device__ void mod_mul_mpz(mpz_t result, const mpz_t a, const mpz_t b, const mpz_t mod)
{
    mpz_mul(result, a, b);
    mpz_mod(result, result, mod);
}

__device__ void gcd_mpz(mpz_t result, const mpz_t a, const mpz_t b)
{
    mpz_gcd(result, a, b);
}

__global__ void pollards_rho_cuda(mpz_t *factor, const mpz_t n, const mpz_t sqrt_n, int *found)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (*found)
        return;

    mpz_t x, y, c, d;
    gmp_randstate_t state;

    gmp_randinit_mt(state);
    gmp_randseed_ui(state, clock64() ^ idx);

    mpz_inits(x, y, c, d, NULL);
    mpz_urandomm(x, state, sqrt_n);
    mpz_add_ui(x, x, 2);
    mpz_set(y, x);
    mpz_urandomm(c, state, sqrt_n);
    mpz_set_ui(d, 1);

    while (mpz_cmp_ui(d, 1) == 0 && !*found)
    {
        mod_mul_mpz(x, x, x, n);
        mpz_add(x, x, c);
        mpz_mod(x, x, n);

        mod_mul_mpz(y, y, y, n);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);
        mod_mul_mpz(y, y, y, n);
        mpz_add(y, y, c);
        mpz_mod(y, y, n);

        mpz_sub(d, x, y);
        mpz_abs(d, d);
        gcd_mpz(d, d, n);
    }

    if (mpz_cmp(d, n) != 0 && mpz_cmp_ui(d, 1) != 0)
    {
        if (atomicExch(found, 1) == 0)
        {
            mpz_set(*factor, d);
        }
    }

    mpz_clears(x, y, c, d, NULL);
    gmp_randclear(state);
}

void pollards_rho_parallel_cuda(mpz_t factor, const mpz_t n, int threads)
{
    mpz_t sqrt_n;
    mpz_init(sqrt_n);
    mpz_sqrt(sqrt_n, n);

    mpz_t *d_factor;
    int *d_found;
    int h_found = 0;

    cudaMalloc(&d_factor, sizeof(mpz_t));
    cudaMalloc(&d_found, sizeof(int));

    cudaMemcpy(d_found, &h_found, sizeof(int), cudaMemcpyHostToDevice);

    int numBlocks = 10;
    pollards_rho_cuda<<<numBlocks, threads>>>(d_factor, n, sqrt_n, d_found);

    cudaDeviceSynchronize();

    cudaMemcpy(&h_found, d_found, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(factor, d_factor, sizeof(mpz_t), cudaMemcpyDeviceToHost);

    cudaFree(d_factor);
    cudaFree(d_found);

    mpz_clear(sqrt_n);
}

int main()
{
    FILE *input_file = fopen("input.txt", "r");
    FILE *output_file = fopen("output_cuda.txt", "w");

    if (!input_file || !output_file)
    {
        printf("Error opening file.\n");
        return 1;
    }

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);

    char n_str[1024];
    int threads = 256;

    while (fscanf(input_file, "%1023s", n_str) != EOF)
    {
        mpz_set_str(n, n_str, 10);

        fprintf(output_file, "Factoring %s using %d CUDA threads...\n", mpz_get_str(NULL, 10, n), threads);

        clock_t start = clock();
        pollards_rho_parallel_cuda(factor, n, threads);
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

    fclose(input_file);
    fclose(output_file);

    mpz_clears(n, factor, NULL);
    return 0;
}
