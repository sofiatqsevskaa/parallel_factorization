#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

long long gcd(long long a, long long b) {
    while (b != 0) {
        long long temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

long long mod_mul(long long a, long long b, long long mod) {
    long long result = 0;
    a = a % mod;
    while (b > 0) {
        if (b % 2 == 1) {
            result = (result + a) % mod;
        }
        a = (a * 2) % mod;
        b /= 2;
    }
    return result;
}

long long pollards_rho(long long n) {
    if (n % 2 == 0) return 2;

    long long x = rand() % n;
    long long y = x;
    long long c = rand() % n;
    long long d = 1;
    int restart = 0; 

    #pragma omp parallel shared(x, y, c, d, restart)
    {
        while (d == 1 && !restart) {
            #pragma omp single
            {
                x = (mod_mul(x, x, n) + c) % n;
                y = (mod_mul(y, y, n) + c) % n;
                y = (mod_mul(y, y, n) + c) % n;
                d = gcd(abs(x - y), n);

                if (d == n) {
                    restart = 1; 
                }
            }
            #pragma omp barrier 
        }
    }

    return d;
}

int main() {
    long long n;

    printf("Enter a composite number: ");
    scanf("%lld", &n);

    if (n <= 1) {
        printf("Invalid input. Enter a number greater than 1.\n");
        return 1;
    }

    long long factor = pollards_rho(n);
    printf("A non-trivial factor of %lld is %lld\n", n, factor);

    return 0;
}
