#include "fp.hpp"
#include "types.hpp"


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

/* Cosas que necesito
 * 		Polinomios en Z, es decir Z[x]
 *		Multiplicar un big_int por un Zx
 *		Comparación del un polinomio en Zx con el número 0
 *		Dividir un polinomio en Zx por un entero
 *		Sumar polinomios en Z_p[x] con polinomios en Z[x]
 *
 * */
bool HenselLifting (const Zx &polynomial, int p, const Fpxelem &u1, const Fpxelem &w1, big_int bound, Zx & u, Zx & w){
	big_int leadCoef = polynomial.lc();
	Zx pol = polynomial * leadCoef;
	Fpxelem::Felem lc(leadCoef);
	u1 *=  (lc * u1.lc().inv()); //This is more efficient than normalize and then multiply by lc
	w1 *=  (lc * w1.lc().inv());

	Fpxelem s, t;
	eea (u1, w1, &s, &t);//This must always be 1. Test it!!

	u = Zx(u1); u[u.getSize()] = leadCoef;
	w = Zx(w1); w[w.getSize()] = leadCoef;
	Zx err = pol - u*w;
	big_int modulus = p;
	bound = 2*bound*leadCoef;

	while (err != 0 && modulus < bound ){
		Fpxelem c(err/modulus);
		pair< Fpxelem, Fpxelem > qr = div2 (s*c, w1);
 		u += Zx(t*c + qr.first * u1) * modulus;
		w += Zx(qr.second) * modulus; 
		err = pol - u*w;
		modulus *= p;
	}

	if (err == 0){
		big_int delta = cont(u);
		u /= delta;
		w /= (leadCoef / delta); //delta must be a divisor of leadCoef (Test it!!)
		return true;
	}
	return false;
}
