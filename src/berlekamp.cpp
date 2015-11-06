#include <vector>
#include <algorithm>
#include "fp.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "types.hpp"


template<typename T>
using matrix = std::vector< std::vector<T> >;


/**
 * Input: a polynomial pol over a field of size q
 * Output: Matrix Q with x^0, x^q, x^{2q},..., x^{(n-1)*q} (mod pol) as rows
 * Complexity: O(q n^2) where n is deg(pol)
 * There is a solution in O(log(q)n^2 + n^3) it is better for big q and small n
 */
template <typename Fxelem>
matrix<typename Fxelem::Felem> formMatrix (const Fxelem &pol) {
    typename Fxelem::F f = pol.getField();
	bint q = f.getSize(), cont = 1;
	int n = pol.deg();
	typename Fxelem::Felem aux = f.get(0);
	std::vector<typename Fxelem::Felem> r(n, f.get(0));
	r[0] = 1; //r == (1, 0, ..., 0)
	matrix<typename Fxelem::Felem> result;
	result.push_back(r);
	for (bint i = 1; i<= (n-1)*q; ++i, ++cont){ //TODO ¿está bien definida la multiplicación (n-1)*q ? (n es un int)
		// r = (-r_{n-1}*pol_0, r_0 -r_{n-1}*pol_1,..., r_{n-2}-r_{n-1}*a_{n-1})
		aux = r[n-1];
		for (ll j = n-1; j >= 1; ++j){
			r[j] = r[j-1]-aux*pol[j];
		}
		r[0] = -aux*pol[0];
		if (cont == q){ //This avoids computing i%q
			result.push_back(r);
			cont = 0;
		}
	}
	return result;
}

/**
 * Input: a square matrix.
 * Output: a basis for the kernel of a matrix. The matrix is destroyed.
 *
 * It forms a lower triangular matrix L. It will satisfies that L^2 = L
 * so as (I-L)L = 0, the non zero rows of I-L form a base for the kernel
 * of the original matrix
 *
 * Complexity:
 *  O(n^3) where n is the dimension of the square matrix
 */
template <typename Fxelem>
std::vector< std::vector< typename Fxelem::Felem > kernelBasis (matrix<typename Fxelem::Felem> & mat){
	bint n = mat.size();
	bint i, j;
	std::vector< std::vector< typename Fxelem::Felem > result;
	for (bint k = 0; k < n; ++k ){
		//Search for pivot element
		for (i = k; i < n && mat[k][i] == 0 ; ++i);

		if (i<n){
			//Normalize column i
			typename Fxelem::Felem inv = mat[k][i].inv();
			for (j = 0; j < n; ++j){
				if (j==k) mat[j][i] = 1; //This is the pivot
				else if (mat[j][i] != 0) mat[j][i] *= inv;
			}
			//Interchange column i with column k
			for (j = 0; j < n; ++j){
				std::swap(mat[j][k], mat[j][i]);
			}
			if (i==k) ++i;
			while(i < n){
				for (j = 0; j < n ; ++j){
					mat[j][i] -= mat[j][k]*mat[k][i];
				}
				++i;
			}
		}
	}
	// M = M - I;
	for (i=0; i<n; ++i)
		mat[i][i] -= 1;
	//Return non zero rows
	j=0;
	while (j < n){
		//Look for the next non zero row
		while (true){
			if (j >=n) break;
			for (i = 0; i < n && mat[j][i] == 0; ++i);
			if (i==n) ++j;
			else break;
		}
		if (j >= n) break;
		result.push_back(mat[j]);//TODO: Quizás se pueda optimizar con un move
	}
	return result; //result[0] should always be (1, 0, ... 0). (Test it!)
}

/* Berlekamp's algorithm
 *
 * Input: a square-free polynomial pol \in F_{p^m}[x]
 * Output: a vector with the irreducible factors of pol
 *
 * Theoretical background:
 *  The set W:={v(x) \in FX | v^q = v (mod pol)} is a vectorial space
 *   whose dimension is the number of irreducible factors of pol. If v \in W
 *   is a non constant polynomial then:
 *    pol(x) = \prod_{s \in F} gcd(v(x)-s, pol(x));
 *  So computing all those gcd where for a base {v_1 .. v_k} of W gives us
 *   the irreducible polynomials of pol
 *
 * Complexity: q is the size of the field and n the degree of pol and k
 * is the number of factors of pol (on average is log(n)):
 *  O(k q n^2 +n^3)
 *
 * */

template <typename Fxelem>
std::vector< Fxelem > berlekamp_simple (const Fxelem &pol){
	std::vector< Fxelem > factors;
	factors.push_back(pol);
	bint r=2;
	auto mat = formMatrix(pol);
	int n = pol.deg();
	for (int i=0; i<n; ++i)
		mat[i][i] -= 1;
//matrix<typename Fxelem::Felem> kernelBasis (const matrix<typename Fxelem::Felem> & mat){
	auto base = kernelBasis<Fxelem>(mat);
	int k = base.size();
	while (factors.size() < k){
		for (int i = 0; i < factors.size(); ++i){
			Fxelem v(base[r]);
			for(auto &s : pol.getField().getElems()){ 
				Fxelem g = gcd(v-s, factors[i]);
				if (g != 1 && g != factors[i]){
					factors[i]/=g; //We continue in the loop with the new factors[i] because it is a divisor of the old factors[i] so it is not necessary to check the previous s and r.
					factors.push_back(g);
					if (factors.size() == k) return factors;
				}
			}
			++r; //TODO: En el libro viene así, no me convence. Mirar.
		}
	}
	return factors;
}

template std::vector< Fpxelem > berlekamp_simple (const Fpxelem &pol);
