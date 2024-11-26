#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*Koristi Evklidov algoritam za generiranje najgolem zaednicki delitel na a i b, taka sto gcd(a,b) = gcd(b,a%b)
Procesot se povtoruva dodeka b ne stane 0, koga b=0, a e najgolemiot zaednicki delitel*/
long long gcd(long long a, long long b) {
    while (b != 0) {
        long long temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

/*Sekvencata koja ja generirame x1...xn se bazira na xi+1 = (xi)^2 + c mod n kaj koj moze da se otkrie ciklus ako 
dve vrednosti vo sekvencata, na primer xi i xj, stanat kongruentni (xi ≡ xj mod p, kade p e faktor na n). 
Otkrivanjeto na takov ciklus ni ovozmozuva da go presmetame GCD(n, |xi - xj|), koj moze da dade netrivijalen delitel na n.*/


/* Funkcijata mod_mul go presmetuva (a * b) mod mod, pritoa ovozmozuva rabota so golemi broevi bez rizik od overflow. */
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



/* 
pollards_rho(n) presmetuva sekvenca koja se generira so formulata: x(i+1) = (x(i)^2 + c) mod n.
Sekvencata pocnuva so sluchajni vrednosti za x i c, i se koristi za da se otkrie ciklus, shto moze da ukaze na nekoi netrivijalni deliteli na n. 
Ako se najde GCD (golem zaednicki delitel) na razlikata od x i y so n, i GCD e različen od 1 i n, togas toj GCD se vrakja kako netrivijalen faktor.
*/
long long pollards_rho(long long n) {
    if (n % 2 == 0) return 2;

    long long x = rand() % n;
    long long y = x;
    long long c = rand() % n;
    long long d = 1;

    while (d == 1) {
        x = (mod_mul(x, x, n) + c) % n;
        y = (mod_mul(y, y, n) + c) % n;
        y = (mod_mul(y, y, n) + c) % n;
        d = gcd(abs(x - y), n);
        if (d == n) return pollards_rho(n);
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
