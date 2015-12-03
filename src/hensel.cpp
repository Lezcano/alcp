#include "fp.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "zelem.hpp"
#include "zxelem.hpp"
#include "types.hpp"
#include "factorizationFq.hpp"
#include "generalPurpose.hpp"
#include <vector>
#include <utility>


//std::vector< pair < Zxelem, unsigned int > > squareFreeFactChar0(const & Zxelem){
//
//}

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
 * 			If there are not such polynomials the function returns false
 */

/* Cosas que necesito
 * 		Polinomios en Z, es decir Z[x]
 *		Multiplicar un big_int por un Zx
 *		Comparación del un polinomio en Zx con el número 0
 *		Dividir un polinomio en Zx por un entero
 *		Sumar polinomios en Z_p[x] con polinomios en Z[x]
 *
 * */
bool HenselLifting (const Zxelem &polynomial, Fpxelem u1, Fpxelem w1, Zxelem & u, Zxelem & w){
	//if (u1.getField.getP() != w1.getField.getP())
	big_int p = u1.getField().getP();
	big_int bound = normInf(polynomial)*fastPow((big_int)2, polynomial.deg());
	big_int leadCoef = polynomial.lc();
	Zxelem pol = polynomial * leadCoef;
	Fpxelem::Felem lc = u1.getField().get(leadCoef);
	u1 *=  (lc * u1.lc().inv()); //This is more efficient than normalize and then multiply by lc
	w1 *=  (lc * w1.lc().inv());

	Fpxelem s(u1.getField().get(0)), t(u1.getField().get(0)); //This is a random value because s & t must be initialized
	eea (u1, w1, s, t);//This must always be 1. Test it!!
	u = Zxelem(u1); u[u.deg()] = leadCoef;
	w = Zxelem(w1); w[w.deg()] = leadCoef;
	Zxelem err = pol - u*w;
	big_int modulus = p;
	bound = 2*bound*leadCoef;

	while (err != 0 && modulus < bound ){
		Fpxelem c(err/modulus, u1.getField().getP());
		auto qr = (s*c).div2(w1);
 		u += Zxelem(t*c + qr.first * u1) * modulus;
		w += Zxelem(qr.second) * modulus;
		err = pol - u*w;
		modulus *= p;
	}

	if (err == 0){
		big_int delta = content(u);
		u /= delta;
		w /= (leadCoef / delta); //delta must be a divisor of leadCoef (Test it!!)
		return true;
	}
	return false;
}
/*
unsigned int heuristic (unsigned int deg, unsigned int numberOfPrimesUsed, const std::vector<unsigned int> & posibilitiesSizes, unsigned int intersectionSize){

}
*/
/*
std::vector< Zxelem > factorizationHenselSquareFree(const Zxelem & poli, HenselSubsets & hs){
	std::vector< Zxelem > result;
	//Preguntar si hs no existe y en tal caso crear una nueva
	HenselSubsets hs();	
	while (hs.oneMorePrime()){
		auto factorsModP = factorizationCantorZassenhaus(Fpxelem (pol, random.getP()));//¿Como coger el primo?
		hs.insert(factorsModP);
	}
	pair < Zxelem, pair < Fpxelem, Fpxelem > > option;
	while (hs.bestOption(option)){	
		Zxelem u(0), w(0);
		if (HenselLifting(option.first, option.second.first, option.second.second, u, w)){
			if (hs.firstIsIrreducible()){
				result.push_back(u);
				hs.removeFirstLastOption();
				if (hs.secondIsIrreducible()){
					result.push_back(w);
					hs.removeSecondLastOption();
				}
				//else we continue in the loop lifting the second
			}	
			else{
				if (hs.secondIsIrreducible()){//we continue in the loop lifting the first
					result.push_back(w);
					hs.removeSecondLastOption();
				}
				else{
					//TODO: copiar hs modificado y llamar a la funcion de nuevo
					result.insert(
						result.end(),
						std::make_move_iterator(aux.begin()),
						std::make_move_iterator(aux.end())
					);
				}
			}	
		}
	}



}

std::vector< pair < Zxelem, unsigned int > factorizationHensel(const Zxelem & pol){
	auto aux = squareFreeFactChar0 (pol);
	std::vector< pair < Zxelem, unsigned int > result;
	for (auto & pair: aux){
		auto factors = factorizationHenselSquareFree (pair.first);
		for (auto & elem : factors){
			result.push_back(elem, pair2.second);
		}
	}
	return result;
}*/
