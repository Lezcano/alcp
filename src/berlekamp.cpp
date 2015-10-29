<<<<<<< HEAD
#include <field.hpp>
#include <vector>
#include <algorithm>

template <typename FX>
typedef std::vector< vector< F > > matrixF;

/**
 * Input: a polinomial pol over a field of size q
 * Output: Matrix Q with x^0, x^q, x^{2q},..., x^{(n-1)*q} (mod pol) as rows
 * Complexity: O(q n^2) where n is deg(pol)
 * There is a solution in O(log(q)n^2 + n^3)
 */
const matrixF&& formMatrix (const FX &pol) {
	bint q = pol.baseFieldSize();
	int n = pol.degree();
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
 * Output: a basis for the kernel of a matrix. The matrix is destroyed.
 */
const vector< vector< F > >&& kernelBasis (const matrixF & matrix){
	bint n = matrix.size();
	bint i, j, k;
	vector< vector< F > > result;
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
	}
	// M = M - I;
	for (i=0; i<n; ++i)
		matrix[i][i] -= 1;
	//Return non zero rows
	j=0;
	while (j < n){
		//Look for the next non zero row
		while (true){
			if (j >=n) break;
			for (i = 0; i < n && matrix[j][i] == 0; ++i);
			if (i==n) ++j;
			else break;
		}	
		if (j >= n) break;
		result.push_back(matrix[j]);
	}
	return result; //result[0] should always be (1, 0, ... 0). (Test it!)
}

//Input: a square-free polynomial pol \in F_{p^m}[x]
const std::vector<FX>&& berlekamp_simple (const FX &pol){
	vector<FX> factors = pol;
	bint r;
	matrixF matrix = formMatrix(pol);
	for (int i=0; i<n; ++i)
		matrix[i][i] -= 1;
	vector< vector< F > > base = kernelBasis(matrix);
	int k = base.size(), j; 
	while (factors.size() < k){
		for (int i = 0; i < factors.size(); ++i){
			for(s \in F){ //Iterar sobre todos los elementos del cuerpo??
				FX g = gcd(FX(base)-s, factors[i]);
				if (g != 1 && g != factors[i]){
					factors[i]/=g; //We continue in the loop with the new factors[i] because it is a divisor of the old factors[i] so it is not necessary to check the previous s and r.
					factors.push_back(g); 
				}
			}
			++r;
		}
	}
}
