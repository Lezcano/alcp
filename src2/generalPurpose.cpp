#include "generalPurpose.hpp"


namespace alcp {

/**
 * Miller Rabin: Primality Test
 *
 * Description:
 *  Montecarlo primality test for positive integers
 *
 * Theoretical background:
 *  A prime number is odd
 *  Z_p is a field iff p is prime
 *  A field does not contain any non-trivial root of the unity
 *  By Fermat's little theorem if p is prime, a^{p-1} = 1 (mod p)
 *  - Decompose n-1 = s*2^r
 *  - For some random a<-{2,..,n-2}
 *  - Check that (a^s)^2^i is not a trivial root of the unity of i=1..r-1 in Z_n
 *      This translates into the test (a^s)^2^i = 1 (mod n)
 *  - Check that (a^s)^2^r = a^{n-1} = 1 in Z_n
 *
 * Complexity:
 *  O(k*log^3(n))
 *  Probability of false positive: 4^{-k}
 *
 * Remarks:
 *  Even though here we implement the classical version of the algorithm, if we
 *   restrict ourselves to long long integers, it is enough to check with:
 *  a = 2,3,4,5,11,13,17,19,23,29,31 and 37
 */
    bool millerRabin(big_int n, int k /*= 35*/) {
        big_int s = n - 1, a;
        int r = 0;

        // Discards edge cases for the random generation process
        if (n == 2 || n == 3) return true;
        // Check that n is odd and greater than 1
        if (n % 2 == 0 || n == 1) return false;

        // Decompose n-1 as n-1=s*2^r
        while (s % 2 == 0) {
            s /= 2;
            r++;
        }
        while (k--) {
            // Generate a random number a<-{2..n-2}
            a = (rand() % (n - 3)) + 2;
            a = fastPowMod(a, s, n);

            // if a=1,n-1, we are not going to find any non trivial root
            // since a^2 (mod n) = 1
            if (a == 1 || a == n - 1) continue;

            for (int i = 0; i < r; i++) {
                a = (a * a) % n;

                // If we find a non trivial root of the unity
                if (a == 1) return false;

                // By fermat's little theorem a^{n-1}=1 (mod n) if n prime
                if (i == r - 1 && a != 1) return false;

                // We are not going to find any nontrivial root of unity
                if (a == n - 1) break;
            }
        }
        return true;
    }



/**
 * Pollard's rho factorization algorithm
 *
 * Explanation:
 *  Given a non-prime number, it returns one of its factors.
 *
 * Theoretical background:
 *  It uses Floyd's cycle detection algorithm to detect
 *   a cycle of the form x_i = y_i (mod p) where p is
 *   the smallest prime number of n.
 *  Floyd's cycle detection algorithm relies on the fact that
 *   if x_{i+1} = f(x_{i}) and x_i does have a cycle, then it
 *   we can find an index j such as x_j = x_{2j}
 *  In the algorithm, x = x_i, y = x_{2i}
 *  The function x \mapsto x^2+1 (mod n) behaves like a random
 *   function for all practical purposes.
 *  The algorithm is based on the fact that both sequences
 *   x_{i+1} = x^2_i + 1(mod n) and x'_{i+1} = x'^2_i + 1 (mod p)
 *   obey the same recurrence and thus, by the birthday paradox,
 *   x_i will find a cycle in O(sqrt(p)) instead of O(sqrt(n)).
 *  Remark: A rigorous anaylsis of Pollard's rho algorithms is
 *   open problem
 *
 * Complexity:
 *  O(sqrt(p)) on average. O(1) bits of space.
 */
    long long pollardRhoBrent(long long n, long long limit /* = 1000 */) {
        long long x = 2, y = 2, p = 1;
        while ((p == 1 || p == n) && limit--) {
            x = (x * x + 1) % n;
            y = (y * y + 1) % n;
            y = (y * y + 1) % n;
            // gcd returns the positive gcd
            p = gcd(y - x, n);
        }
        return p;
    }

    // Suppose n prime, else r^-1 should be computed more carefully
    bool pollardRhoLogarithm(long long g, long long h, long long n, long long &log) {
        n--;
        auto f = [=](long long &x, long long &a, long long &b) {
            switch (x % 3) {
                   case 0: x = x*x % (n+1);  a =  a*2  % n;  b =  b*2  % n;  break;
                   case 1: x = x*g % (n+1);  a = (a+1) % n;                  break;
                   case 2: x = x*h % (n+1);                  b = (b+1) % n;  break;
            }
        };
        long long x1 = 1, a1 = 0, b1 = 0;
        long long x2 = 1, a2 = 0, b2 = 0;
        // Floyd's cycle algorithm to find a collision
        // Invariant g^a1*h^b1 = g^a2*h^b2
        do {
            f(x1, a1, b1);
            f(x2, a2, b2);
            f(x2, a2, b2);
        } while (x2 != x1);
        long long r = (b2 - b1) % n;
        if(r < 0)
            r += n;
        if (r == 0) return false;
        long long ir, aux, aa;
          aa = (a1 - a2) % n;
        if(aa < 0)
            aa += n;
        long long gc = gcd(r, aa);
        r /= gc; aa /= gc;
        if(eea(r, n, ir, aux) != 1)
            return false;
        ir %= n;
        if(ir < 0)
            ir += n;
        log = (ir * aa) % n;
        return true;
    }


}
