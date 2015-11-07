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
template <typename Fpxelem, Zx>
/* Cosas que necesito
 * 		Polinomios en Z, es decir Z[x]
 *		Multiplicar un bint por un Zx
 *		Comparación del un polinomio en Zx con el número 0
 *		Dividir un polinomio en Zx por un entero
 *		Sumar polinomios en Z_p[x] con polinomios en Z[x]
 *
 * */
/*
 * Nota: pol lo voy a modificar, pero luego lo voy a dejar como estaba
 *
 * */
pair<Fpxelem, Fpxelem> HenselLifting (Zx &pol, int p, const Fpxelem &u1, const Fpxelem &w1, bint bound){//Lo suyo sería devolver un struct...
	bint leadCoef = pol.lc();
	pol *= leadCoef;
	Fpxelem::Felem lc(leadCoef);
	u1 *=  (lc * u1.lc().inv()); //This is more efficient than normalize and then multiply by lc
	w1 *=  (lc * w1.lc().inv());

	Fpxelem s, t;
	eea (u1, w1, &s, &t);//This must always be 1. Test it!!

	Zx u(u1); u[u.getSize()] = leadCoef;
	Zx w(w1); w[w.getSize()] = leadCoef;
	Zx err = pol - u*w;
	bint modulus = p;
	bound = 2*bound*leadCoef;
	while (err != 0 && modulus < bound ){
		Fpxelem c(err/modulus);
		pair< Fpxelem, Fpxelem > qr = div2 (s*c, w1);
 		u += (t*c + qr.first * u1)*modulus; //Quiero que la multiplicación me dé un polinomio en Z[x]
		w += qr.second * modulus; //Tambien necesito que el resultado de la multiplicación esté en Z[x]
		err = pol - u*w;
		modulus *= p;
	}

	if (err == 0){
		bint delta = cont(u);
		u /= delta;
		w /= (leadCoef / delta); //delta must be a divisor of leadCoef
		return pair<u, w>;
	}
	return pair<NULL, NULL>;
}
