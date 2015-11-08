#include <vector>
#include "fp.hpp"
#include "types.hpp"

/* Auxiliary function. Finds x such that
 * ax = 1 (mod q)
 * */
bint reciprocal (bint a, bint q){
	bint x, y;
	eea (a, q, x, y);
	return x;
}

/**
 * Integer Chinese Remainder Algorithm (Garner's Algorithm)
 * Description:
 *  Given positive moduli m_i \in Z (0 \leq i \leq n) which are
 *   relatively prime and given corresponding residues u_i \in Z_{m_i}
 *   compute the unique integer u \in Z_m (where m = \prod m_i) such that
 *   u = u_i (mod m_i) i = 0,...,n
 *
 * Theoretical background:
 *  The algorithm express the solution in the mixed radix representation:
 * 		u = v_0 + v_1 m_0 + ... + v_n (\prod_{i=0}^{n-1} m_i)
 *  where v_k \in Z_{m_k} k = 0..n
 *
 *  Now v_0 = u_0 and for k >= 1 we have (mod m_k):
 * 		u_k = v_0 + ... + v_k (\prod_{i=0}^{k-1})
 * 	and we know everything but v_k so we compute it solving the equation
 * 	(we need to compute the inverse of \prod_{i=0}^{k-1} (mod m_k). We
 * 	precompute all those inverses.
 *
 * Complexity:
 *  O(n^2)
 *
 */
const bint& integerCRA (const std::vector<bint> & m, const std::vector<bint> & u){
	int n = m.size()-1;
	bint prod, aux;
	std::vector<bint> inv, v;
	for (int k = 1; k <= n; ++k ){
		for	(int i = 0 ; i <= k-1; ++i)
			prod = (prod*m[i])%m[k];
		inv.push_back(reciprocal(prod, m[k]));	//inv starts in 0, in the book it starts in 1
	}
	v.push_back(u[0]);
	for (int k = 1; k <= n; ++k){
		aux = v[k-1];
		for (int j = k-2; j >=0; --j)
			aux = (aux*m[j]+v[j])%m[k];
		v[k] = ((u[k]-aux)*inv[k-1])%m[k];	//it is inv[k-1] because inv starts in 0
	}
	bint result = v[n];
	for (int k = n-1; k >= 0; --k ){
		result = result*m[k] + v[k];
	}
	return result;
}
