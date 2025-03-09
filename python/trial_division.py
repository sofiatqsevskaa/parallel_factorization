import numpy as np
import time
from numba import cuda
import math


@cuda.jit
def trial_division_kernel(n, factors):
    idx = cuda.grid(1)
    if idx == 0:
        if n % 2 == 0:
            factors[0] = 2
            return

    factor = 3 + 2 * idx
    if factor * factor > n:
        return

    if n % factor == 0:
        factors[0] = factor


def trial_division_gpu(n, max_factor=None):
    if n % 2 == 0:
        return 2

    if max_factor is None:
        max_factor = int(math.sqrt(n)) + 1

    num_tests = (max_factor - 3) // 2 + 1
    factors = np.zeros(1, dtype=np.int64)
    d_factors = cuda.to_device(factors)

    threads_per_block = 256
    blocks_per_grid = (num_tests + threads_per_block - 1) // threads_per_block

    trial_division_kernel[blocks_per_grid, threads_per_block](n, d_factors)
    factors = d_factors.copy_to_host()

    return factors[0] if factors[0] > 0 else None


def factorize_from_file(input_filename, output_filename):
    with open(input_filename, "r") as infile, open(output_filename, "w") as outfile:
        overall_start = time.time()

        for line in infile:
            n = int(line.strip())
            outfile.write(f"Factoring {n}...\n")
            outfile.flush()

            start = time.time()
            factor = trial_division_gpu(n)
            end = time.time()

            print(f"Attempting to factor {n}...")

            if factor:
                outfile.write(f"Non-trivial factor: {factor}\n")
                print(f"Result for {n}: {factor}")
            else:
                outfile.write("Failed to find a non-trivial factor.\n")
                print("No factor")

            outfile.write(f"Time taken: {end - start:.6f} seconds\n\n")
            outfile.flush()

        overall_end = time.time()
        outfile.write(f"Overall runtime: {
                      overall_end - overall_start:.6f} seconds\n\n")


factorize_from_file("input.txt", "output_trial_division.txt")
