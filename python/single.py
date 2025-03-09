import random
import time
from sympy import gcd, sqrt


def mod_mul(a, b, mod):
    return (a * b) % mod


def pollards_rho(n):
    if n % 2 == 0:
        return 2

    sqrt_n = int(sqrt(n))

    x = random.randrange(1, sqrt_n + 1)
    y = x
    c = random.randrange(1, sqrt_n + 1)
    d = 1

    while d == 1:
        x = (mod_mul(x, x, n) + c) % n
        y = (mod_mul(y, y, n) + c) % n
        y = (mod_mul(y, y, n) + c) % n
        d = gcd(abs(x - y), n)

    if d is None or d == n or d == 1:
        return None

    return d if d <= sqrt_n else n // d  # Ensure factor ≤ sqrt(n)


def factorize_from_file(input_filename, output_filename):
    with open(input_filename, "r") as infile, open(output_filename, "w") as outfile:
        overall_start = time.time()

        for line in infile:
            n = int(line.strip())
            outfile.write(f"Factoring {n}...\n")

            start = time.time()
            factor = pollards_rho(n)
            end = time.time()

            print(f"Attempting to factor {n}...")

            if factor and factor != n and factor != 1:
                outfile.write(f"Non-trivial factor: {factor}\n")
                print(f"Non-trivial factor: {factor}\n")
            else:
                outfile.write("Failed to find a non-trivial factor.\n")
                print("Failed to find a non-trivial factor.\n")

            outfile.write(f"Time taken: {end - start:.6f} seconds\n\n")

        overall_end = time.time()
        outfile.seek(0)
        outfile.write(f"Overall runtime: {
                      overall_end - overall_start:.6f} seconds\n\n")


factorize_from_file("input.txt", "golemiot.txt")
