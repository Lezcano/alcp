#include "fp.hpp"
#include "types.hpp"


#define Felem typename Fxelem::Felem
#define F typename Fxelem::F

/**
 * Input:
 *		1) A primitive polynomial pol
 *		2) A prime integer p which does not divide pol[pol.degree]
 *		3) Two relatively primes polynomials u1, w1 \in Z_p[X] such that
 *			pol = u1*w1 (mod p)
 *		4) A bound B of all integer coefficients of pol and all coefficients
 *			of u, v
 *
 * Output: Two polynomials u, w (if there are such polynomials) such that
 * 			pol = u*w and such that u' = u1' (mod p) and w' = w1' (mod p)
 * 			where v' stands for the monic normalization of v (mod p)
 * 			If there are not such polynomials ...
 */
template <typename Fxelemi, Zx>
/* Cosas que necesito
 * 		Polinomios en Z, es decir Z[x]
 *		Multiplicar un int por un Zx
 *		En las asignaciones de u1 y w1 de las primeras lineas hago dos cosas equivalentes, la primera sintacticamente correcta, estaría bien poder hacer la segunda
 *		Comparación del un polinomio en Zx con el número 0
 *		Dividir un polinomio en Zx por un entero
 *
 * */
pair<Fxelem, Fxelem> HenselLifting (const Zx &pol, int p, const Fxelem &u1, const Fxelem &w1, bint bound){//Lo suyo sería devolver un struct...
	ll leadCoef = ;
	Felem lc(pol.lc());
	pol *= pol.lc();
	u1 *=  (lc * u1.lc().inv());
	w1 *=  (lc * w1.lc().inv());

	Fxelem s, t;
	eea (u1, w1, &s, &t);//This must always be 1. Test it!!

	Zx u(u1); u[u.getSize()] = leadCoef;
	Zx w(w1); w[w.getSize()] = leadCoef;
	Zx err = pol - u*w;
	bint modulus = p;
	bound = 2*bound*leadCoef;
	while (err != 0 && modulus < bound ){
		Fxelem c(err/modulus);

		pair< Fxelem, Fxelem > qr = div2 (s*c, w1);


	}

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
const vector< vector< Felem > >&& kernelBasis (const matrix & mat){
	bint n = mat.size();
	bint i, j;
	vector< vector< Felem > > result;
	for (bint k = 0; k < n; ++k ){
		//Search for pivot element
		for (i = k; i < n && mat[k][i] == 0 ; ++i);

		if (i<n){
			//Normalize column i
			Felem inv = mat[k][i].inv();
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
		result.push_back(mat[j]);
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
 *  whose dimension is the number of irreducible factors of pol. If v \in W
 *  is a non constant polynomial then:
 *   pol(x) = \prod_{s \in F} gcd(v(x)-s, pol(x));
 *  So computing all those gcd where for a base {v_1 .. v_k} of W gives us
 *  the irreducible polynomials of pol
 *
 * Complexity: q is the size of the field and n the degree of pol and k
 * is the number of factors of pol (on average is log(n)):
 *  O(k q n^2 +n^3)
 *
 * */
const std::vector< Fxelem >&& berlekamp_simple (const Fxelem &pol){
	vector< Fxelem > factors = pol;
	bint r;
	matrix mat = formMatrix(pol);
	for (int i=0; i<n; ++i)
		mat[i][i] -= 1;
	vector< vector< Felem > > base = kernelBasis(mat);
	int k = base.size(), j;
	while (factors.size() < k){
		for (int i = 0; i < factors.size(); ++i){
			for(auto &s : F.getElems()){ //Iterar sobre todos los elementos del cuerpo??
				Fxelem g = gcd(Fxelem(base)-s, factors[i]);
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
