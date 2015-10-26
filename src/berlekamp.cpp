#include <field.hpp>
#include <vector>

template <typename FX>

//Input: a polinomial pol over a field of size q
//O(q n^2) where n is deg(pol)
const matrixF& formMatrix (const FX &pol) { //TODO: change/check ALL ll
	ll q = pol.baseFieldSize(), n = pol.degree(); 
	FX::F aux(0, q);
	std::vector<F> r(n, 0);
	r[0] = 1; //r == (1, 0, ..)
	matrixF result;
	result.push_back(r);
	for (ll i = 1; i<= (n-1)*q; ++i){
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

//Input: a square-free polynomial pol \in F_{p^m}[x]
const std::vector<FX>& berlekamp_simple (const FX &pol){

}
