#include <cstdlib> // abs function for long long ints
#include "exceptions.hpp"
#include "types.hpp"
#include "fpxelem.hpp"

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
template<typename T>
T fastPowMod(T a, ll b, T p){
    if(b==0)return 1;
    if(b%2 != 0){
    	return (a*fastPowMod(a,b-1,p))%p;
    }
    else{
    	ll aux = fastPowMod(a,b/2,p);
    	return (aux*aux)%p;
    }
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
bool millerRabin(ll n, int k /*= 35*/){
    ll s=n-1, r=0, a;

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
 * Extended Euclidean Algorithm
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
ll eea(ll a, ll b, ll& x, ll& y){
    // We set a = |a| and b = |b| and set a flag if the sign was changed
    bool ca = false, cb = false;
    if(a < 0)a = -a, ca = true;
    if(b < 0)b = -b, cb = true;

    x = 1; y = 0;
    ll xx = 0, yy = 1;
    ll q, r, r1, r2;
    while(b){
        // The following invariant holds:
        //  a = x*|a|+y*|b|
        //  b = xx*|a|+yy*|b|

        // Compute quotient and reminder
        q = a/b; r = a-q*b;

        // After this r = r1*|a|+r2*|b| holds
        r1 = x-q*xx; r2 = y-q*yy;

        // Iterate
        a = b; x = xx; y = yy;
        b = r; xx = r1; yy = r2;
    }
    if(ca) x = -x;
    if(cb) y = -y;

    return std::abs(a);
}

ll gcd(ll a, ll b){
    ll x,y;
    return eea(a,b,x,y);
}

/**
 * Extended Euclidean Algorithm for an arbitrary DFU
 */
template<typename T>
T eea (T a, T b, T &x, T &y){
    typename T::F f = a.getField();
    if(a == 0){
        if(b == 0)
            throw EOperationUnsupported("Cannot compute the greatest common divisor of two zero elements.");
        x = f.get(0);
        y = unit(a);
        return normalForm(b);
    }
    if(b == 0){
        x = unit(a);
        y = f.get(0);
        return normalForm(a);
    }
    a = normalForm(a);
    b = normalForm(b);

    T ua = unit(a), ub = unit(b);

      x = T(f.get(1)); y = T(f.get(0));
    T xx = T(f.get(0)), yy = T(f.get(1));
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

template Fpxelem eea (Fpxelem a, Fpxelem b, Fpxelem &x, Fpxelem &y);

template<typename T>
T gcd(T a, T b){
    typename T::F f = a.getField();
    T x = f.get(0), y = f.get(0);
    return eea(a,b,x,y);
}
template Fpxelem gcd(Fpxelem a, Fpxelem b);
