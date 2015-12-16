#include "hensel.hpp"

#include <vector>
#include <utility>
#include <iostream> // TODO Quitar

#include "fp.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "zelem.hpp"
#include "zxelem.hpp"
#include "types.hpp"
#include "henselSubsets.hpp"
#include "factorizationFq.hpp"
#include "generalPurpose.hpp"

const bool verbose = false;

//std::vector< pair < Zxelem_b, unsigned int > > squareFreeFactChar0(const & Zxelem_b){
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
bool HenselLifting (const Zxelem_b &polynomial, Fpxelem_b u1, Fpxelem_b w1, Zxelem_b & u, Zxelem_b & w){
	//TODO: if (u1.getField.getP() != w1.getField.getP())
	big_int p = u1.getField().getP();
	big_int bound = normInf(polynomial)*fastPow((big_int)2, polynomial.deg());
	big_int leadCoef = polynomial.lc();
	Zxelem_b pol = polynomial * leadCoef;
	Fpxelem_b::Felem lc = u1.getField().get(leadCoef);
	u1 *=  (lc * u1.lc().inv()); //This is more efficient than normalize and then multiply by lc
	w1 *=  (lc * w1.lc().inv());

	Fpxelem_b s(u1.getField().get(0)), t(u1.getField().get(0)); //This is a random value because s & t must be initialized
	eea (u1, w1, s, t);//This must always be 1. Test it!!
	u = Zxelem_b(u1); u[u.deg()] = leadCoef;
	w = Zxelem_b(w1); w[w.deg()] = leadCoef;
	Zxelem_b err = pol - u*w;
	big_int modulus = p;
	bound = 2*bound*leadCoef;

	while (err != 0 && modulus < bound ){
		Fpxelem_b c(err/modulus, u1.getField().getP());
		auto qr = (s*c).div2(w1);
 		u += Zxelem_b(t*c + qr.first * u1) * modulus;
		w += Zxelem_b(qr.second) * modulus;
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

std::vector< Zxelem_b > factorizationHenselSquareFree(Zxelem_b poli, HenselSubsets & hs){
	std::vector< Zxelem_b > result;
	int asd =0, primes[13] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43};
	while (hs.oneMorePrime()){
		/*/
		 big_int prime = randomPrime();//TODO Coger el codigo del gcd modular y ponerlo con una cota más grande¿?
		/*/
		if (verbose && asd >= 13) {
			std::cout << "Los primos no valen para este polinomio" << std::endl;
			return result;
		}
		big_int prime = primes[asd++];
		/**/
		Fpxelem_b aux(poli, prime);
		if (poli.lc()%prime == 0 || gcd(aux, Fpxelem_b(poli.derivative(), prime)) != 1)
			continue;

		auto factorsModP = factorizationCantorZassenhaus(aux);
		hs.insert(factorsModP, aux);
	}

	Option option = hs.bestOption();
	while (option.b){
		Zxelem_b u(0), w(0);
        if(verbose){
            std::cout << "Intento elevar esto:" << std::endl;
            std::cout << option.u << std::endl << option.w << std::endl;
        }
		if (HenselLifting(poli, option.u, option.w, u, w)){
			//The first factor is always irreducible because hs first iterates through the options in ascending order with respect to the degree of the first polynomial
			result.push_back(u);
			poli = poli/u;
            if(verbose){
                std::cout << "¡Ha funcionado!" << std::endl;
                std::cout << "u:" <<  u <<std::endl;
                std::cout << "w:" <<  w <<std::endl;
            }
			hs.removeFirstLastOption(w);
		}
        if(verbose){
            std::cout << "=====" << std::endl;
        }
		option = hs.bestOption();
	}
	Zxelem_b last = hs.getLast();
	if (last != 1)
		result.push_back(last);
	return result;
}

std::vector< Zxelem_b > factorizationHenselSquareFree(const Zxelem_b & poli){
	HenselSubsets hs(poli);
	return factorizationHenselSquareFree(poli, hs);
}
/*
std::vector< std::pair < Zxelem_b, unsigned int > > factorizationHensel(const Zxelem_b & pol){
	//auto aux = squareFreeFactChar0 (pol);
	std::vector< std::pair < Zxelem_b, unsigned int > result;
	for (auto & pair: aux){
		auto factors = factorizationHenselSquareFree (pair.first);
		for (auto & elem : factors){
			result.push_back(elem, pair.second);
		}
	}
	return result;
}
*/
