#include <field.hpp>
#include <vector>
#include <algorithm>

template <typename FX>

/**
 * Input: a polinomial pol over a field of size q
 * O(q n^2) where n is deg(pol)
 */
const matrixF&& formMatrix (const FX &pol) { //TODO: change/check ALL ll
	bint q = pol.baseFieldSize();
	int n = pol.degree(); //TODO: n es el grado de un polinomio, ¿usar int?
	FX::F aux(0, q);
	std::vector<F> r(n, 0);
	r[0] = 1; //r == (1, 0, ..., 0)
	matrixF result;
	result.push_back(r);
	for (bint i = 1; i<= (n-1)*q; ++i){ //TODO ¿está bien definida la multiplicación (n-1)*q ? (n es un int)
		// r = (-r_{n-1}*pol_0, r_0 -r_{n-1}*pol_1,..., r_{n-2}-r_{n-1}*a_{n-1})
		aux = r[n-1];
		for (ll j = n-1; j>= 1; ++j){
			r[j] = r[j-1]-aux*pol[j];
		}
		r[0] = -aux*pol[0];
		if (i%q == 0)
			result.push_back(r);
	}	
	return result;	
}
/**
 * Input: a square matrix.
 * Output: a basis for the kernel of matrix. The matrix is destroyed.
 */
const vector< vector< F > >&& kernelBasis (const matrixF & matrix){
	bint n = matrix.size();
	bint i, j, k;
	for (k = 0; k < n; ++k ){
		//Search for pivot element
		for (i = k; i<n && matrix[k][i] ==0 ; ++i);

		if (i<n){
			//Normalize column i
			F inv = matrix[k][i].inv();
			for (j = 0; j < n; ++j){
				if (j==k) matrix[j][i] = 1; //This is the pivot
				else if (matrix[j][i] != 0) matrix[j][i]*= inv;
			}
			//Interchange column i with column k
			F* aux;
			for (j = 0; j < n; ++j){
				std::swap(matrix[j][k], matrix[j][i]);
			}
			if (i==k) ++i;
			while(i < n){
				for (j = 0; j < n ; ++j){
					matrix[j][i] -= matrix[j][k]*matrix[k][i]; 
				}	
				++i;
			}
		}
		// M = M - I;
		for (i=0; i<n; ++i){
			matrix[i][i] -= 1;
		}
		//Return non zero rows
		j=0;
		while (j < n){
			while (true){
				for (){
					
				}
			}	
		}
	}
}

//Input: a square-free polynomial pol \in F_{p^m}[x]
const std::vector<FX>& berlekamp_simple (const FX &pol){
	
}
