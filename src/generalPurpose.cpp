#include <cstdlib> // abs function for long long ints
#include "generalPurpose.hpp"
#include "exceptions.hpp"
#include "types.hpp"
#include "zelem.hpp" // getZero y tal
#include "fpxelem.hpp"
#include "fqxelem.hpp"
#include "fpelem.hpp"
#include "fqelem.hpp"
#include "fp.hpp"
#include "fq.hpp"

#include <iostream> //TODO quitar

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
T fastPowMod(const T&a, U b, const T& p){
    T aux = a;
    T result = getOne(a);
    while (b != 0){
        if (b % 2 == 0){
            aux = (aux*aux)%p;
            b /= 2;
        }
        else{
            result = (result*aux)%p;
            b -= 1;
        }
    }
    return result;
}

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
 *      This tranlsates into the test (a^s)^2^i = 1 (mod n)
 *  - Check that (a^s)^2^r = a^{n-1} = 1 in Z_n
 *
 * Complexity:
 *  O(k*log^3(n))
 *  Probability of false posiive: 4^{-k}
 *
 * Remarks:
 *  Even though here we implement the classical version of the algorithm, if we
 *   restrict ourselves to long long integers, it is enough to check with:
 *  a = 2,3,4,5,11,13,17,19,23,29,31 and 37
 */
bool millerRabin(big_int n, int k /*= 35*/){
    big_int s=n-1, a;
    int r=0;

    // Discards edge cases for the random generation process
    if(n==2 || n==3) return true;
    // Check that n is odd and greater than 1
    if(n%2==0 || n==1) return false;

    // Decompose n-1 as n-1=s*2^r
    while(s%2==0){
        s/=2;
        r++;
    }
    while(k--) {
        // Generate a random number a<-{2..n-2}
        a=(rand()%(n-3))+2;
        a=fastPowMod(a, s, n);

        // if a=1,n-1, we are not going to find any non trivial root
        // since a^2 (mod n) = 1
        if(a==1 || a == n-1) continue;

        for(int i=0; i<r; i++){
            a=(a*a)%n;

            // If we find a non trivial root of the unity
            if(a==1) return false;

            // By fermat's little theorem a^{n-1}=1 (mod n) if n prime
            if(i==r-1 && a!=1) return false;

            // We are not going to find any nontrivial root of unity
            if(a==n-1) break;
        }
    }
    return true;
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
T eea (T a, T b, T &x, T &y){
    if(a == 0){
        if(b == 0)
            throw EOperationUnsupported("Cannot compute the greatest common divisor of two zero elements.");
        x = getZero(a);
        y = unit(b);
        return normalForm(b);
    }
    if(b == 0){
        x = unit(a);
        y = getZero(b);
        return normalForm(a);
    }
    a = normalForm(a);
    b = normalForm(b);
    T ua = unit(a), ub = unit(b);

      x = getOne(a); y = getZero(b);
    T xx = getZero(a), yy = getOne(b);
    while(b!=0){
        // The following invariant holds:
        //  a = x*|a|+y*|b|
        //  b = xx*|a|+yy*|b|
        // Compute quotient and reminder
        T q = a/b, r = a-q*b;

        // After this r = r1*|a|+r2*|b| holds
        T r1 = x-q*xx, r2 = y-q*yy;

        // Iterate
        a = b; x = xx; y = yy;
        b = r; xx = r1; yy = r2;
    }
    x /= unit(ua*unit(a));
    y /= unit(ub*unit(a));
    return normalForm(a);
}


template<typename T>
T gcd(T a, T b){
    T aux = getZero(a);
    if (a == 0 || b == 0 ) return getOne(a);
    while(b != 0){
        aux = b;
        b = a%b;
        a = aux;
    }
    return normalForm(a);
}


template<typename T, typename U>
T fastPow (const T& a, U b){
    T aux = a;
    T result = getOne(a);
    while (b != 0){
        if (b % 2 == 0){
            aux*=aux;
            b /= 2;
        }
        else{
            result *= aux;
            b -= 1;
        }
    }
    return result;
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
long long pollardRhoBrent (long long n){
    long long x=2, y=2, p=1;
    while(p==1){
        x = (x*x+1)%n;
        y = (y*y+1)%n;
        y = (y*y+1)%n;
        // gcd returns the positive gcd
        p = gcd(y-x, n);
    }
    return p;
}

// Suppose n prime, else the r^-1 should be computed more carefully
bool pollardRhoLogarithm(long long g, long long h, long long n, long long & log){
    long long x1 = 1, a1 = 0, b1 = 0;
    long long x2 = 1, a2 = 0, b2 = 0;
    auto f = [h, g, n](long long &x, long long &a, long long &b){
        switch(x%3){
            case 0: x*=h;            b = (b+1)%n;  break;
            case 1: x*=x; a=(a*2)%n; b = (b*2)%n;  break;
            case 2: x*=g; a=(a+1)%n;               break;
        }};
    // Floyd's cycle algorithm to find a collision
    // Invariant g^a1*h^b1 = g^a2*h^b2
    do{
        f(x1,a1,b1);
        f(x2,a2,b2);
        f(x2,a2,b2);
    }while(x2!=x1);
    long long r = (b1 - b2) % n;
    if(r==0) return false;
    long long ir, aux;
    // compute r1 = r^-1
    eea(r,n,ir,aux);
    log = (ir*(a2-a1))%n;
    return true;
}



template big_int fastPowMod<big_int, big_int>(const big_int &a, big_int b, const big_int & p);
template Fpxelem fastPowMod<Fpxelem, big_int>(const Fpxelem &a, big_int b, const Fpxelem & p);

template big_int fastPow<big_int, int>(const big_int& a, int b);
template big_int fastPow<big_int, unsigned int>(const big_int& a, unsigned int b);
template big_int fastPow<big_int, unsigned long long>(const big_int& a, unsigned long long b);
//template big_int fastPow<big_int, std::size_t>(const big_int& a, std::size_t b);
template big_int fastPow<big_int, big_int>(const big_int& a, big_int b);
template Fpelem fastPow<Fpelem, big_int>(const Fpelem& a, big_int b);
template Fqelem fastPow<Fqelem, big_int>(const Fqelem& a, big_int b);

template big_int gcd<big_int>(big_int a, big_int b);
template Fpxelem gcd<Fpxelem>(Fpxelem a, Fpxelem b);
template Fqxelem gcd<Fqxelem>(Fqxelem a, Fqxelem b);

template big_int eea<big_int>(big_int a, big_int b, big_int &x, big_int &y);
template Fpxelem eea<Fpxelem>(Fpxelem a, Fpxelem b, Fpxelem &x, Fpxelem &y);
template Fqxelem eea<Fqxelem>(Fqxelem a, Fqxelem b, Fqxelem &x, Fqxelem &y);
