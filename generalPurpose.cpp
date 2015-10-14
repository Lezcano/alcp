
#include "Exceptions.h"

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
ll fastPowMod(ll a, ll b, ll p){
    if(b==0)return 1;
    if(b%2){
        ll aux = fastPowMod(a,b/2,p);
        return (b*b)%p;
    }
    return (a*fastPowMod(a,b-1,p))%p;
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
bool millerRabin(ll n, int k=35) {
	ll s=n-1, r=0, a, aux;

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
        // Generate a random number a<-{1..n-2}
		a=(rand()%(n-2))+1;
		a=fastPowMod(a, s, n);

        // if a=1,n-1, we are not going to find any non trivial root
        // since a^2 (mod n) = 1
        if(a==1 || a == n-1) continue;

		for(int i=0; i<r; i++){
            a=(a*a)%n;

            // If we find a non trivial root of the unity
            if(a==1) return false;

            // We are not going to find any nontrivial root of unity
			if(a==n-1) break;

            // By fermat's little theorem a^{n-1}=1 (mod n) if n prime
            if(i==r-1 && a!=1) return false;

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
 *   It returns d
 *
 * Theoretical background:
 *  The algorithm is based in the equality mcd(a,b)=mcd(b%a,a)
 *  At the end of each recursive call we have the invariant
 *   a*x+b*y=gdc(a,b)
 *
 * Complexity:
 *  O(log(min(a,b)))
 *
 */
ll eea (ll a, ll b, ll& x, ll& y) {
    // If a=0 => gcd(0,b) = b = 0*0+1*b
    if(a==0){
        x=0; y=1;
        return b;
    }
    ll x1, y1;

    // Recursively compute the gcd
    ll gcd = eea(b%a, a, x1, y1);

    // By induction hypothesis gcd = b%a*x1+a*y1
    // So we update the coefficients so that gcd = a*x+b*y
    x= y1 - (b/a)*x1;
    y = x1;

    return gcd;
}


#endif // __GENERAL_PURPOSE_H
