#ifndef __GENERAL_PURPOSE_H
#define __GENERAL_PURPOSE_H

#include "types.hpp"
// idk why this can't be on the .cpp
#include "zelem.hpp"

namespace alcp {
    /**
     * Exponentiation by Squaring
     *
     * Explanation:
     *  Computes a^b (mod m)
     * Theoretical background:
     *  The algorithm is based on the fact that:
     *  a^0 = 1
     *  a^{2k} (mod m) = (a^k (mod m) * a^k (mod m)) (mod m)
     *  a^{2k+1} (mod m) = (a*(a^{2k} (mod m)) (mod m)
     *
     * Complexity:
     *  O((b*log(a))) supposing multiplication in O(1)
     */
    template<typename T, typename U>
    T fastPowMod(const T &a, U b, const T &p) {
        T aux = a;
        T result = getOne(a);
        while (b != 0) {
            if (b % 2 == 0) {
                aux = (aux * aux) % p;
                b /= 2;
            }
            else {
                result = (result * aux) % p;
                b -= 1;
            }
        }
        return result;
    }

    /**
     * Extended Euclidean Algorithm for an arbitrary DE
     *
     * Description:
     *  Given two integers a, b it computes:
     *   ax+by=d=gcd(a,b)
     *   It returns d with  d > 0, i.e. in its normal form.
     *
     * Theoretical background:
     *  The algorithm is based in the equality mcd(a,b)=mcd(b%a,a)
     *  The details of the invariants are commented in the code
     *  This algorithm works with a few tweaks for an arbitrary Euclidean Domain
     *
     * Complexity:
     *  O(log(min(a,b)))
     *
     */
    template<typename T>
    T eea(T a, T b, T &x, T &y) {
        T zero = getZero(a);
        T one = getOne(a);
        if (a == 0) {
            if (b == 0) {
                x = zero;
                y = zero;
                return zero;
            }
            x = zero;
            y = unit(b);
            return normalForm(b);
        }
        if (b == 0) {
            x = unit(a);
            y = zero;
            return normalForm(a);
        }
        a = normalForm(a);
        b = normalForm(b);
        T ua = unit(a), ub = unit(b);

        x = one;
        y = zero;
        T xx = zero, yy = one;
        while (b != 0) {
            // The following invariant holds:
            //  a = x*|a|+y*|b|
            //  b = xx*|a|+yy*|b|
            // Compute quotient and reminder
            T q = a / b, r = a - q * b;

            // After this r = r1*|a|+r2*|b| holds
            T r1 = x - q * xx, r2 = y - q * yy;

            // Iterate
            a = b;
            x = xx;
            y = yy;
            b = r;
            xx = r1;
            yy = r2;
        }
        x /= unit(ua * unit(a));
        y /= unit(ub * unit(a));
        return normalForm(a);
    }

    template<typename T>
    T gcd(T a, T b) {
        T aux;
        if (a == 0 && b == 0)
            return a;
        while (b != 0) {
            aux = b;
            b = a % b;
            a = aux;
        }
        return normalForm(a);
    }


    template<typename T, typename U>
    T fastPow(const T &a, U b) {
        T aux = a;
        T result = getOne(a);
        while (b != 0) {
            if (b % 2 == 0) {
                aux *= aux;
                b /= 2;
            }
            else {
                result *= aux;
                b -= 1;
            }
        }
        return result;
    }

    bool millerRabin(big_int n, int k = 35);

    long long pollardRhoBrent(long long n, long long limit = 1000);

    bool pollardRhoLogarithm(long long g, long long h, long long n, long long &log);
}
#endif // __GENERAL_PURPOSE_H
