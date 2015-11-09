#include <vector>
#include <algorithm>
#include <utility>
#include "fp.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "generalPurpose.hpp"
#include "types.hpp"
#include <iostream> //TODO Quitar

/* Detalles de la implementación:
 * 	formMatrix:
 * 		-En la función  el for interno se recorre en sentido descendente para
 * 		poder hacer los calculos sin usar otro array.
 * 		-Se lleva un contador en el for externo que cuenta hasta q para no
 * 		tener que usar %
 * 	kernelBasis:
 * 		-El for de la linea 86 va en sentido descendente porque se necesita en
 * 		todo momento el valor mat[k][i]. De esta manera éste se actualiza
 * 		exactamente en la ultima operación.
 * 		-El for de la línea 94 debería calcular I-M pero eso es muy	caro y
 * 		nosotros sólo necesitamos una base, así que en lugar de eso, calculo 
 * 		M-I que también vale como base.
 * 		- Linea 98 j=1 se hace (en vez de j=0) porque el primer elemento
 * 		de la base es siempre (1, 0,..,0), (aunque como yo cojo la
 * 		base con los números opuestos sería (-1, 0,..,0)) y no se usa
 * 		para nada así que directamente no la calculo ni añado a result
 * 		(esto hace que en berlekamp r se inicialize a 0 en vez de
 * 		a 1 y que k se inicialize a base.size()+1
 * 		-Línea 103, la fila j es nula si y solo si el elemento mat[j][j]
 * 		es cero (en caso contrario es 1)
 *  berlekamp_simple:
 *  	-Lo dicho antes sobre las inicializaciones de r y k
 *	partialFactorDD:
 *		-Resulta que para elevar (en mod pol) un polinomio a la q, lo
 *		unico que hay que hacer es multiplicar sus coeficientes por
 *		la matriz de formMatrix, así que para calcular x^{iq}-x (mod pol) 
 *		lo que hago es cogerme el x^{(i-1)q} que tenía de antes, multiplico
 *		sus	coeficientes por la matriz y ya tengo x^{iq} (mod pol)
 *		-La primera iteración no la hago dentro del bucle porque x^q es la
 *		segunda fila de la matriz, así que no tengo que calcularlo.
 *		-En el result.push_back de dentro del while hago 
 *		gcd(x^{qi}-x(mod pol), pol1) donde pol es el polinomio original y pol1
 *		es un divisor. En el libro hacen modulo pol1 en vez de pol, pero notese que al ser pol multiplo de pol1 se tiene que
 *		(x^{qi}-x(mod pol)) (mod pol1) = x^{qi}-x(mod pol1)
 *		así que por la propiedad del algorimo de euclides el gcd es el mismo.
 *
 *		*/

template<typename T>
using matrix = std::vector< std::vector<T> >;


/**
 * Input: a polynomial pol over a field of size q
 * Output: Matrix Q with x^0, x^q, x^{2q},..., x^{(n-1)*q} (mod pol) as rows
 * Complexity: O(q n^2) where n is deg(pol)
 * There is a solution in O(log(q)n^2 + n^3), it is better for big q and small n
 */
template <typename Fxelem>
matrix<typename Fxelem::Felem> formMatrix (const Fxelem &pol) {
    typename Fxelem::F f = pol.getField();
	bint q = f.getSize(), cont = 1;
	int n = pol.deg();
	//typename Fxelem::Felem aux = f.get(0);
	std::vector<typename Fxelem::Felem> r(n, f.get(0));
	r[0] = 1; //r == (1, 0, ..., 0)
	matrix<typename Fxelem::Felem> result;
	result.push_back(r);
	for (bint i = 1; i<= (n-1)*q; ++i, ++cont){ //TODO ¿está bien definida la multiplicación (n-1)*q ? (n es un int)
		// r = (-r_{n-1}*pol_0, r_0 -r_{n-1}*pol_1,..., r_{n-2}-r_{n-1}*a_{n-1})
		auto aux = r[n-1];
		for (ll j = n-1; j >= 1; --j){
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
std::vector< std::vector< typename Fxelem::Felem > > kernelBasis (matrix<typename Fxelem::Felem> & mat){
	bint n = mat.size();
	bint i, j;
	std::vector< std::vector< typename Fxelem::Felem > > result;

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
			if (i!=k){
				for (j = 0; j < n; ++j){
					std::swap(mat[j][k], mat[j][i]);
				}
			}
			i=0;
			while(i < n){
				if (i==k) {
					++i;
					continue;
				}
				for (j = n-1; j >= k ; --j){//It has to be backwards
					mat[j][i] -= mat[j][k]*mat[k][i];
				}
				++i;
			}
		}
	}
	// M = M - I; //Note this is -1*(I-M)
	for (i=0; i<n; ++i)
		mat[i][i] -= 1;

	//Return non zero rows
	j=1; //we do not need the first row
	while (j < n){
		//Look for the next non zero row
		while (true){
			if (j >=n) break;
			if (mat[j][j] == 0) ++j; //The row is zero iff mat[j][j] == 0
			else break;
		}
		if (j >= n) break;
		result.push_back(mat[j]);//TODO: Quizás se pueda optimizar con un move
		++j;
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
	bint r=0;
	auto mat = formMatrix(pol);
	int n = pol.deg();
	for (int i=0; i<n; ++i)
		mat[i][i] -= 1;
	auto base = kernelBasis<Fxelem>(mat);
	int k = base.size()+1;//we do not have computed the first element of the base, so have to add 1 to k
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

/*
//Part I
template<typename Fxelem>
std::vector< Fxelem > squareFreeFactorization (const Fxelem &pol);
 */

//Part II
/*
 *
 *
 * */
template<typename Fxelem>
std::vector< std::pair< Fxelem, unsigned int> > partialFactorDD ( Fxelem &pol){//TODO: Copiar el polinomio en vez de pasarlo por referencia??
	int n = pol.deg();
	auto mat = formMatrix(pol);

	//first iteration is performed out of the loop because we have r in mat (there is no need to compute it again)
	std::vector<typename Fxelem::Felem> r = mat[1];
	//result[i].first will be a product of irreducible polynomials with degree result[i].second
	std::vector< std::pair< Fxelem, unsigned int > > result;
	unsigned int i = 1;
	r[1] -= 1;
	result.push_back(std::make_pair(gcd(Fxelem(r), pol), i));

	r[1] += 1;
	if (result.back().first != 1)
		pol /= result.back().first;

	++i;
	while (i <= pol.deg()/2){
		std::vector<typename Fxelem::Felem> aux = r;
		for (int j = 0; j < n; ++j){
				r[j] = aux[0]*mat[0][j] ;
			for (int k = 1; k < n; ++k ){
				r[j] += aux[k]*mat[k][j];
			}
		}//This is just r = r*mat;
		r[1] -= 1;
		result.push_back(std::make_pair(gcd(Fxelem(r), pol), i+1));//gcd (a_1, w (mod a)) = gcd (a_1, w (mod a_1)) where a_1 divides a (because (w (mod a))(mod a_1) = w (mod a_1))
		r[1] += 1;
		if (result.back().first != 1)
			pol /= result.back().first;
		else
			result.pop_back();
		++i;
	}
	if (pol != 1)
		result.push_back(std::make_pair(pol, pol.deg()));

	return result;	
}
/*
//Part III
template<typename Fxelem>
std::vector< Fxelem > splitFactorsDD (const Fxelem &pol){
}


template std::vector< Fxelem > splitFactorsDD (const Fxelem &pol);
template std::vector< Fxelem > partialFactorDD (const Fxelem &pol);

*/
template std::vector< std::pair< Fpxelem, unsigned int> > partialFactorDD ( Fpxelem &pol);
template std::vector< Fpxelem > berlekamp_simple (const Fpxelem &pol);
