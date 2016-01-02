#include <vector>
#include "types.hpp"
#include "generalPurpose.hpp"

namespace alcp {
/* Auxiliary function. Finds x such that
 * ax = 1 (mod q)
 * */
	big_int reciprocal(big_int a, big_int q) {
		big_int x, y;
		eea(a, q, x, y);
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
	big_int integerCRA(const std::vector<big_int> &m, const std::vector<big_int> &u) {
		const int n = static_cast<int>(m.size() - 1);
		big_int prod, aux;
		std::vector<big_int> inv(n), v(n + 1);
		if (m.size() != u.size())
			throw std::runtime_error("The vectors in integerCRA have to be of the same size.");
		for (int k = 1; k <= n; ++k) {
			prod = m[0] % m[k];
			for (int i = 1; i <= k - 1; ++i)
				prod = (prod * m[i]) % m[k];
			inv[k - 1] = reciprocal(prod, m[k]);
		}
		v[0] = u[0];
		for (int k = 1; k <= n; ++k) {
			aux = v[k - 1];
			for (int j = k - 2; j >= 0; --j)
				aux = (aux * m[j] + v[j]) % m[k];
			v[k] = ((u[k] - aux) * inv[k - 1]) % m[k];
			// Convert to symmetric representation
			if (v[k] > m[k] / 2)
				v[k] -= m[k];
			else if (v[k] < -m[k] / 2)
				v[k] += m[k];
		}
		// Compute the result using horner's algorithm
		big_int result = v[n];
		for (int k = n - 1; k >= 0; --k) {//This must be int because size_t is always >=0 so the loop would be infinite
			result = result * m[k] + v[k];
		}
		return result;
	}
}
