#include "bchCodes.hpp"

#include <vector>
#include <utility>
#include <random>
#include <chrono>
#include <set>
#include <string>
#include <iostream>
#include <unordered_map>

#include "berlekampMassey.hpp"
#include "generalPurpose.hpp"
#include "types.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "fqelem.hpp"
#include "fqxelem.hpp"

namespace alcp {
	std::pair<std::set<int>, Fpxelem_b > BCH::randomErrors(Fpxelem_b v){
        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> distr(0, maxErrors);
        std::mt19937 generator2(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> distr2(0, v.deg());
        std::mt19937 generator3(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> distr3(0, static_cast<int>(g.getField().getP()-1));

        int actualNumErrors = distr(generator);
        std::set<int> indices_set;
		if (v == getZero(v[0]))
			return std::make_pair(indices_set, v);

        for (int i = 0; i < actualNumErrors; i++){
        	int aux, aux2;
        	do{
        		aux = distr2(generator2);
        	}while(indices_set.find(aux) != indices_set.end());//Do not repeat index
        	do{
				aux2 = distr3(generator3);
			}while(v[aux] == aux2);//An error must be a real error
        	indices_set.insert(aux);
            v[aux] = v[0].getField().get(aux2);
        }
        return std::make_pair(indices_set, v);
	}

	/* It checks if the multiplicative order of q (mod l) is m
	 * */
	bool checkOrderMod(const big_int & q, size_t l, size_t m){
		big_int aux = 1;
		for(unsigned int i = 0 ; i < m-1 ; i++){
			aux = (aux*q)%l;
			if (aux == 1)
				return false;
		}
		aux = (aux*q)%l; // here: aux == q^m
		if (aux == 1)
			return true;
		else
			return false;
	}
	/* This solves the linear system Ax = b for a non singular square matrix with no zeros in the main diagonal
	 */
	std::vector<Fqelem_b> solve(std::vector<std::vector<Fqelem_b> > &m, std::vector<Fqelem_b> & b) {
		int i, j, k, rows = m.size(), cols = rows;
		Fqelem_b c;
		for (i = 0; i < rows; i++) {
			b[i] /= m[i][i];
			for (j = cols-1; j >= 0; j--)
				m[i][j] /= m[i][i];
			for (j = i+1; j < rows; j++){
				b[j] -= b[i] * m[j][i];
				for(k = cols-1; k >= i; --k){
					m[j][k] -= m[i][k] * m[j][i];
				}
			}
		}//Gauss finished
		//Now solve triangular system with 1's in the main diagonal (now it is not necessary to edit the matrix)
		for (j = cols-1; j >= 0; --j){
			for(i = j-1; i >= 0 ; --i)
				b[i] -= b[j] * m[i][j];
		}
		return b;
	}

	Fpxelem_b generating_polynomial(const Fqelem_b & alpha, size_t c, size_t d, const big_int & q){
		std::set<Fqelem_b>  rootSet;
		Fqelem_b root = fastPow(alpha, c);
		for (size_t i = 0; i <= d-2; i++){//d is always >= 2
			if (rootSet.find(root) != rootSet.end()){ //Continue if the root was already processed
				root *= alpha;
				continue;
			}
			Fqelem_b aux = root;
			do{
				rootSet.insert(aux);
				aux = fastPow(aux, q ); //Frobenius automorphism
			}while(root != aux);
			root *= alpha;
		}
		std::vector<Fqelem_b> monomial(2, getZero(alpha));
		monomial[0] = getOne(alpha);
		Fqxelem_b result(monomial);
		monomial[1] = getOne(alpha);
		for(auto & r : rootSet){
			monomial[0] = -r;
			result *= Fqxelem_b(monomial);
		}
		//The result's coefficients are in Fp
		std::vector<Fpelem_b> vec(result.deg()+1);
		for(unsigned int i = 0; i <= result.deg(); i++ ){
			vec[i] = static_cast<Fpxelem_b>(result[i])[0]; //coger el elemento de Fp de este elemento de Fq
		}
		return Fpxelem_b(vec);
	}

	/* Precondition: primitive_poly must be a primitive polynomial over a finite field with p elements, with p prime.
	 * The primitivity is not checked.
	 *
	 * An object to encode and decode a BCH code is created. The generating polynomial will be over a finite field of
	 * p^n elements, where p is the prime size of the field such that primitive_poly is over and n is a positive integer.
	 *
	 * The code will be in F_{p^n}[X]/< x^n - 1>.
	 * Let alpha be a primitive l-root of the unity that lives in an extension F_{p^{mn}}.
	 * m is computed as the multiplicative order of q:= p^n module l. Thus, we need to check that gcd(l, p)==1
	 * The primitive_poly provided must have order m*n and it will be used to generate the extension F_{p^{mn}}.
	 *
	 * The generating polynomial of the code will be the lcm of the minimum polynomials of alpha^c, alpha^{c+1},
	 *  \dots, alpha^{c+d-2}.
	 * We can compute the roots of the minimum polynomial of alpha^i applying iteratively the Frobenius automorphism
	 * \phi: x \mapsto x^q
	 *
	 * Note that without any generalization we always have n = 1. Right now the method is implemented for n = 1.
	 * */
	BCH::BCH(const Fpxelem_b &primitive_poly, size_t n, size_t l, size_t cc, size_t d):
		field_ext(primitive_poly){
		//if (primitive_poly.deg() % n != 0)
		//	throw EBadBCHInitialization("The degree of the primitive polynomial is not coherent.");
		this->c = cc;
		big_int p = primitive_poly[0].getField().getP();
		big_int q = fastPow(p, n); //fastPow is unnecessary if n==1
		size_t m = primitive_poly.deg()/n;
		if (! ( 2 <= d && d <= l && gcd(big_int(l), p) == 1 && checkOrderMod(q, l, m) ))
			throw EBadBCHInitialization("Bad Initialization. Check that 2 <= d <= l and that the degree m of the primitive polynomial is the multiplicative order of p (mod l).");

		std::vector<Fpelem_b> aux(2, getZero(primitive_poly[0])); //gen
		aux[1] = 1;
		Fpxelem_b aux2 = Fpxelem_b(aux);
		alpha = field_ext.get(fastPow(aux2, (fastPow(q, m)-1)/l )); // alpha := x^{(q^m-1)/l} \in F_{p^{mn}}
		g = generating_polynomial(alpha, c, d, q); //Esto en el ordenador va a estar como un polinomio sobre F_{q^m} pero
														//sus coeficientes van a estar en realidad sobre F_q (q podría ser p en este momento)
		//g always divides x^l-1, we don't have to take module x^l-1
		dimension = l - g.deg()-1;
		distance = d;
		maxErrors = (d-1)/2;
		length = l;
		std::cout << "generating polynomial: " << g << std::endl << std::endl;
	}

	Fpxelem_b BCH::encode(Fpxelem_b message){
		if(message.deg()+1 > dimension)
			throw EBadFormatMessageBCH("The message is too long.");

		return g*message;
	}

	Fpxelem_b BCH::decode(Fpxelem_b w){
		std::cout << std::endl << std::endl << "Decoding. Computing syndromes." << std::endl;
		std::vector<Fqelem_b> syndromes(distance-1);
		Fqelem_b aux = fastPow(alpha, c);
		//inmersión de w en el anillo de polinomios de la extension F_q
		std::vector<Fpxelem_b> inter(w.deg()+1);
		for( size_t i = 0; i<= w.deg(); i++){
			inter[i] = Fpxelem_b(w[i]);
		}
		Fqxelem_b ww = Fqxelem_b(inter, field_ext); //This is the natural inmersion of w \in F_p to ww \in F_q
		//Fqxelem_b ww = toFqxelem(inter, field_ext); //This is the natural inmersion of w \in F_p to ww \in F_q
		for (size_t i = 0; i <= distance-2; i++ ){
			syndromes[i] = ww.eval(aux);
			aux *= alpha;
		}
		std::cout << "Decoding using Berlekamp algorithm."<< std::endl;
		Fqxelem_b errorLocatorPoly = berlekampMassey<Fqxelem_b>(syndromes);
		std::cout << "Berlekamp finished, computing the roots indices."<< std::endl;
		//std::cout << "errorLocatorPoly " << errorLocatorPoly << std::endl;
		unsigned int nErrors = errorLocatorPoly.deg();
		if(nErrors > maxErrors)
			throw ETooManyErrorsBCH("Too many errors");

		size_t i = 0, index = 1;
		Fqelem_b aux2 = alpha;
		std::vector<int> pos_errors(nErrors);
		while (i != nErrors && index <= length){
			if (errorLocatorPoly.eval(aux2) == 0){
				pos_errors[i++] = length - index; //This is because we work with the reciprocal polynomial
			}
			index++;
			aux2 *= alpha;
		}
		std::cout << nErrors <<  " error/s has/have been detected at index/indices ";
		for (auto &elem : pos_errors){
			std::cout << elem << ", ";
		}
		std::cout << std::endl;
		std::vector<std::vector<Fqelem_b> > m(nErrors, std::vector<Fqelem_b>(nErrors));
		std::vector<Fqelem_b> b(nErrors);
		for (size_t j = 0; j < nErrors; j++){
			Fqelem_b alpha_i = fastPow(alpha, pos_errors[j]); //alpha^{i_1}
			m[0][j] = fastPow(alpha, c*pos_errors[j]); //alpha^{c*i_j}
			for (size_t i = 1; i < nErrors; i++){
				m[i][j] = m[i-1][j]*alpha_i;
			}
		}

		for (size_t i = 0; i < nErrors; i++)
			b[i] = syndromes[i];

/*
		for (size_t i = 0; i < nErrors; i++){
			for (size_t j = 0; j < nErrors; j++){
				std::cout << m[i][j] << "   ";
			}
			std::cout << b[i] << std::endl;
		}
*/

		solve(m, b); //The solution is in b. m is changed

		//Now we proceed to correct the errors
		for (size_t i = 0; i < nErrors; i++){
			w[pos_errors[i]] -= static_cast<Fpxelem_b>(b[i])[0]; //Although b[i] is a Fqelem, it is actually a Fpelem
		}
		return w;
	}
	Fpxelem_b BCH::getG() const { return g; }

	size_t BCH::getDimension() const { return dimension; }
}
