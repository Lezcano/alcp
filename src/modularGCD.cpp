//EN PROCESO... NO LO LEAS XD

big_int intContent (const FX & a){
	//Calcula el gcd de los coeficientes con el signo del coeficiente director
}


/* Modular GCD 
 *  Given A, B \in Z[x_1..x_k] nonzero, it obtains gcd (A, B) via modular
 *  reduction.
 * */
//TODO: ¿Como representamos los polynomios en varias variables?
const modularGCD FX & (const FX & a, const FX & b){
	big_int ia = intContent(a); a /= ia;
	big_int ib = intContent(b); b /= ib;
	//Compute coefficient bound of gcd(a, b)
	big_int c = gcd (ia, ib);
	big_int g = gcd (a.leadingCoef, b.leadingCoef); //TODO
	if (g < 0) g = -g;
	big_int q = 0, h = 0, p;
	int n = min (a.degree(), b.degree());
	big_int limit = (1<<n)*g*min(normInf(a), normInf(b));
	
	while (true){
		p = generateNewLargePrime();
		while (p %g == 0)
			p = generateNewLargePrime();

	}
}
