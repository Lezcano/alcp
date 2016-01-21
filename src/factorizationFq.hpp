#ifndef __FACTORIZATION_FQ
#define __FACTORIZATION_FQ

#include <vector>
#include <utility>
#include <algorithm>
#include <utility>
#include <random>
#include <chrono>

#include "fpxelem.hpp"
#include "fqxelem.hpp"
#include "generalPurpose.hpp"
#include "types.hpp"
namespace alcp {
    /* Detalles de la implementación:
     *  formMatrix:
     *      -En la función  el for interno se recorre en sentido descendente para
     *      poder hacer los calculos sin usar otro array.
     *      -Se lleva un contador en el for externo que cuenta hasta q para no
     *      tener que usar %
     *  kernelBasis:
     *      -El for de la linea 86 va en sentido descendente porque se necesita en
     *      todo momento el valor mat[k][i]. De esta manera éste se actualiza
     *      exactamente en la ultima operación.
     *      -El for de la línea 94 debería calcular I-M pero eso es muy caro y
     *      nosotros sólo necesitamos una base, así que en lugar de eso, calculo
     *      M-I que también vale como base.
     *      - Linea 98 j=1 se hace (en vez de j=0) porque el primer elemento
     *      de la base es siempre (1, 0,..,0), (aunque como yo cojo la
     *      base con los números opuestos sería (-1, 0,..,0)) y no se usa
     *      para nada así que directamente no la calculo ni añado a result
     *      (esto hace que en berlekamp r se inicialize a 0 en vez de
     *      a 1 y que k se inicialize a base.size()+1
     *      -Línea 103, la fila j es nula si y solo si el elemento mat[j][j]
     *      es cero (en caso contrario es 1)
     *  berlekamp_simple:
     *      -Lo dicho antes sobre las inicializaciones de r y k
     *  partialFactorDD:
     *      -Resulta que para elevar (en mod pol) un polinomio a la q, lo
     *      unico que hay que hacer es multiplicar sus coeficientes por
     *      la matriz de formMatrix, así que para calcular x^{iq}-x (mod pol)
     *      lo que hago es cogerme el x^{(i-1)q} que tenía de antes, multiplico
     *      sus coeficientes por la matriz y ya tengo x^{iq} (mod pol)
     *      -La primera iteración no la hago dentro del bucle porque x^q es la
     *      segunda fila de la matriz, así que no tengo que calcularlo.
     *      -En el result.push_back de dentro del while hago
     *      gcd(x^{qi}-x(mod pol), pol1) donde pol es el polinomio original y pol1
     *      es un divisor. En el libro hacen modulo pol1 en vez de pol, pero notese que al ser pol multiplo de pol1 se tiene que
     *      (x^{qi}-x(mod pol)) (mod pol1) = x^{qi}-x(mod pol1)
     *      así que por la propiedad del algorimo de euclides el gcd es el mismo.
     *
     *      */

    template<typename T>
    using matrix = std::vector<std::vector<T> >;

    template<typename Fxelem>
    matrix<typename Fxelem::Felem> formMatrix(const Fxelem &pol);

//Part I
//Outputs a vector of pairs with the factors and multiplicities of the square free factorization (not necesarily sorted by multiplicity)
/*
 *
 * */
    template<typename Fxelem>
    std::vector<std::pair<Fxelem, unsigned int> > squareFreeFF(Fxelem a) {
        unsigned int i = 1;
        std::vector<std::pair<Fxelem, unsigned int> > result;

        Fxelem b = a.derivative();
        if (b != 0) {
            Fxelem c = gcd(a, b);
            Fxelem w = a / c;
            while (w != 1) {
                Fxelem y = gcd(w, c);
                Fxelem z = w / y;
                if (z.deg() != 0 )
                    result.push_back(std::make_pair(z, i));
                i++;
                w = y;
                c = c / y;
            }
            if (c != 1) {
                big_int p = c.getField().getP();
                // c is now of the form: c_0 + c_p x^p + ... c_{kp} x^{kp}
                big_int exponent = fastPow(p, c.getField().getM() - 1);
                //This for computes c = c^{1/p}
                std::vector<typename Fxelem::Felem> rootPOfC;
                for (std::size_t j = 0; j <= c.deg(); j += static_cast<std::size_t>(p))
                    if (c[j] != 0)
                        rootPOfC.push_back(fastPow(c[j], exponent));
                    else {
                        rootPOfC.push_back(getZero(c[j]));
                    }
                auto aux = squareFreeFF(Fxelem(rootPOfC));

                for (auto &pair: aux) {
                    pair.second *= static_cast<unsigned int>(p);
                }
                result.insert(
                        result.end(),
                        std::make_move_iterator(aux.begin()),
                        std::make_move_iterator(aux.end())
                );
            }
        }
        else { // a is of the form: a_0 + a_p x^p + ... a_{kp} x^{kp} (because its derivative is zero)
            big_int p = a.getField().getP();
            big_int exponent = fastPow(p, a.getField().getM() - 1);
            //This for computes a = a^{1/p}
            std::vector<typename Fxelem::Felem> rootPOfA;
            for (int j = 0; j <= a.deg(); j += static_cast<unsigned int>(p))
                if (a[j] != 0)
                    rootPOfA.push_back(fastPow(a[j], exponent));
                else {
                    rootPOfA.push_back(getZero(a[j]));
                }
            auto aux = squareFreeFF(Fxelem(rootPOfA));

            for (auto &pair: aux) {
                pair.second *= static_cast<unsigned int>(p);
            }
            result.insert(
                    result.end(),
                    std::make_move_iterator(aux.begin()),
                    std::make_move_iterator(aux.end())
            );
        }
        return result;
    }

    //Part II
    /*
     *
     *
     * */
    template<typename Fxelem>
    std::vector<std::pair<Fxelem, unsigned int> > partialFactorDD(Fxelem pol) {
    	//result[i].first will be a product of irreducible polynomials with degree result[i].second
    	std::vector<std::pair<Fxelem, unsigned int> > result;

    	int n = pol.deg();
    	if (n == 1){
    		result.push_back(std::make_pair(pol, 1));
    		return result;
    	}
        //auto mat = formMatrix(pol);
		auto mat = formMatrixBigQ(pol);

        //first iteration is performed out of the loop because we have r in mat (there is no need to compute it again)
        std::vector<typename Fxelem::Felem> r = mat[1];


        unsigned int i = 1;
        r[1] -= 1;
        result.push_back(std::make_pair(gcd(Fxelem(r), pol), i));
        r[1] += 1;
        if (result.back().first != 1)
            pol /= result.back().first;
        else
            result.pop_back();

        ++i;
        while (i <= pol.deg() / 2) {
            std::vector<typename Fxelem::Felem> aux = r;
            for (int j = 0; j < n; ++j) {
                r[j] = aux[0] * mat[0][j];
                for (int k = 1; k < n; ++k) {
                    r[j] += aux[k] * mat[k][j];
                }
            }//This is just r = r*mat;
            r[1] -= 1;
            result.push_back(std::make_pair(gcd(Fxelem(r), pol),
                                            i));//gcd (a_1, w (mod a)) = gcd (a_1, w (mod a_1)) where a_1 divides a (because (w (mod a))(mod a_1) = w (mod a_1))
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

    template<typename Fxelem>
    void fastPowModPol(Fxelem &a, big_int b, std::vector<Fxelem> pwrsX, int deg) {
        if (b == 0) {
            a = getOne(a);
        }
        else {
            Fxelem aux = a;
            a = getOne(a);
            while (b != 0) {
                if (b % 2 == 0) {
                    aux *= aux;
                    for (int i = 0; i <= (int) (aux.deg()) - deg; ++i) {//aux.deg is always <= 2*deg-2
                        if (aux[i + deg] != 0)
                            aux += Fxelem(aux[i + deg]) * pwrsX[i];
                    }
                    b /= 2;
                }
                else {
                    a *= aux;
                    for (int i = 0; i <= (int) (a.deg()) - deg; ++i) {//a.deg is always <= 2*deg-2
                        if (a[i + deg] != 0)
                            a += Fxelem(a[i + deg]) * pwrsX[i];
                    }
                    b -= 1;
                }
            }
        }
    }


    template<typename Fxelem>
    Fxelem randomPol(const typename Fxelem::F &field, int degree) {
        //TODO esto genera números aleatorios de 64, parece suficiente, pero si q es mayor que 2^63 en realidad no lo es...
        std::vector<typename Fxelem::Felem> r;
        std::mt19937_64 generator(std::chrono::system_clock::now().time_since_epoch().count());
        for (int i = 0; i <= degree; ++i) {
            r.push_back(field.get(generator()));
        }
        return Fxelem(r);
    }

//Part III
/*
 *
 * */
    template<typename Fxelem>
    std::vector<Fxelem> splitFactorsDD(Fxelem pol, int n) {
        int polDeg = pol.deg();
        if (polDeg <= n) {
            std::vector<Fxelem> factors;
            factors.push_back(pol);
            return factors;
        }
        int m = polDeg / n;

        std::vector<Fxelem> pwrsX;
        std::vector<typename Fxelem::Felem> r(2 * polDeg - 1, getZero(pol.lc()));
        r[polDeg - 1] = 1; //r == (0, 0, ..., 1)
        for (int i = polDeg; i <= 2 * polDeg - 2; ++i) {
            // r = (-r_{n-1}*pol_0, r_0 -r_{n-1}*pol_1,..., r_{n-2}-r_{n-1}*pol_{n-1})
            auto aux = r[polDeg - 1];
            for (std::size_t j = polDeg - 1; j >= 1; --j) {
                r[j] = r[j - 1] - aux * pol[j];
            }
            r[0] = -aux * pol[0];
            r[i] = -1;
            pwrsX.push_back(Fxelem(r));
            r[i] = 0;
        }

        while (true) {
            Fxelem v = randomPol<Fxelem>(pol.getField(), 2 * n - 1);
            if (pol.getField().getSize() % 2 == 0) {//size %2 == 0 iff p %2 == 0
                Fxelem aux = v;
                for (int i = 1; i <= n * m - 1; ++i) {
                    aux *= aux;
                    //This loop performs the operation (mod pol)
                    for (int i = polDeg; i <= aux.deg(); ++i) {//aux.deg is always <= 2*polDeg-2
                        if (aux[i] != 0)
                            aux += Fxelem(aux[i]) * pwrsX[i - polDeg];
                    }
                    v += aux;
                }
            }
            else {
                fastPowModPol<Fxelem>(v, (fastPow(pol.getField().getSize(), n) - 1) / 2, pwrsX, pol.deg());
                v -= getOne(pol);
            }
            Fxelem g = gcd(pol, v);
            if (g != 1 && g != pol) {
                std::vector<Fxelem> factors = splitFactorsDD(g, n);
                std::vector<Fxelem> factors2 = splitFactorsDD(pol / g, n);
                factors.insert(
                        factors.end(),
                        std::make_move_iterator(factors2.begin()),
                        std::make_move_iterator(factors2.end())
                );
                return factors;
            }
        }
    }


/**
 * Input: a polynomial pol over a field of size q
 * Output: Matrix Q with x^0, x^q, x^{2q},..., x^{(n-1)*q} (mod pol) as rows
 * Complexity: O(q n^2) where n is deg(pol)
 * There is a solution in O(log(q)n^2 + n^3), it is better for big q and small n
 */
    template<typename Fxelem>
    matrix<typename Fxelem::Felem> formMatrix(const Fxelem &pol) {
        big_int q = pol.getField().getSize();
        int n = pol.deg();

        std::vector<typename Fxelem::Felem> r(n, getZero(pol.lc()));
        r[0] = 1; //r == (1, 0, ..., 0)
        matrix<typename Fxelem::Felem> result;
       
        result.push_back(r);
        for (big_int i = 1; i <= (n - 1) * q; ++i) { 
            // r = (-r_{n-1}*pol_0, r_0 -r_{n-1}*pol_1,..., r_{n-2}-r_{n-1}*pol_{n-1})
            auto aux = r[n - 1];
            for (std::size_t j = n - 1; j >= 1; --j) {
                r[j] = r[j - 1] - aux * pol[j];
            }
            r[0] = -aux * pol[0];
            if (i % q == 0)
                result.push_back(r);
        }
        return result;
    }
	template<typename Fxelem>
    matrix<typename Fxelem::Felem> formMatrixBigQ(const Fxelem &pol) {
		int polDeg = pol.deg();
		std::vector<typename Fxelem::Felem> rr(polDeg, getZero(pol.lc()));
        rr[0] = 1; //rr == (1, 0, ..., 0)
		matrix<typename Fxelem::Felem> result;
        result.push_back(rr);
		if (polDeg == 1) 
			return result;

		std::vector<Fxelem> pwrsX;

        std::vector<typename Fxelem::Felem> r(2 * polDeg - 1, getZero(pol.lc()));
        r[polDeg - 1] = 1; //r == (0, 0, ..., 1)
        for (int i = polDeg; i <= 2 * polDeg - 2; ++i) {
            // r = (-r_{n-1}*pol_0, r_0 -r_{n-1}*pol_1,..., r_{n-2}-r_{n-1}*pol_{n-1})
            auto aux = r[polDeg - 1];
            for (int j = polDeg - 1; j >= 1; --j) {
                r[j] = r[j - 1] - aux * pol[j];
            }
            r[0] = -aux * pol[0];
            r[i] = -1;
            pwrsX.push_back(Fxelem(r));
            r[i] = 0;
        }
       Fxelem xq = Fxelem({getZero(pol.lc()), getOne(pol.lc())});
       fastPowModPol<Fxelem>(xq, pol.getField().getSize(), pwrsX, polDeg);
       Fxelem aux = xq;
       auto aux2 = static_cast<std::vector<typename Fxelem::Felem> >(xq);
       aux2.resize(polDeg, getZero(pol.lc()));
       result.push_back(aux2); //x^q mod pol

        for (int i = 2; i <= (int)(polDeg) - 1; ++i) {
			 aux = aux*xq; //(x^{i*q});
			for (int j = 0; j <= (int) (aux.deg()) - polDeg; ++j) {//At this point result[i].deg is always <= 2*deg-2
       			if (aux[j + polDeg] != 0)
       	      		aux += Fxelem(aux[j + polDeg]) * pwrsX[j];
      		 }
			auto aux2 = static_cast<std::vector<typename Fxelem::Felem> >(aux);
		    aux2.resize(polDeg, getZero(pol.lc()));
		    result.push_back(aux2);
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
    template<typename Fxelem>
    std::vector<std::vector<typename Fxelem::Felem> > kernelBasis(matrix<typename Fxelem::Felem> &mat) {
        int n = mat.size();
        int i, j;
        std::vector<std::vector<typename Fxelem::Felem> > result;

        for (int k = 0; k < n; ++k) {
            //Search for pivot element
            for (i = k; i < n && mat[(size_t) k][i] == 0; ++i);

            if (i < n) {
                //Normalize column i
                typename Fxelem::Felem inv = mat[k][i].inv();
                for (j = 0; j < n; ++j) {
                    if (j == k) mat[j][i] = 1; //This is the pivot
                    else if (mat[j][i] != 0) mat[j][i] *= inv;
                }
                //Interchange column i with column k
                if (i != k) {
                    for (j = 0; j < n; ++j) {
                        std::swap(mat[j][k], mat[j][i]);
                    }
                }
                i = 0;
                while (i < n) {
                    if (i == k) {
                        ++i;
                        continue;
                    }
                    for (j = n - 1; j >= k; --j) {//It has to be backwards
                        mat[j][i] -= mat[j][k] * mat[k][i];
                    }
                    ++i;
                }
            }
        }
        // M = M - I; //Note this is -1*(I-M)
        for (i = 0; i < n; ++i)
            mat[i][i] -= 1;

        //Return non zero rows
        j = 1; //we do not need the first row
        while (j < n) {
            //Look for the next non zero row
            while (true) {
                if (j >= n) break;
                if (mat[j][j] == 0) ++j; //The row is zero iff mat[j][j] == 0
                else break;
            }
            if (j >= n) break;
            result.push_back(mat[j]);
            ++j;
        }
        return result;

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
    template<typename Fxelem>
    std::vector<Fxelem> berlekamp_simple(const Fxelem &pol) {
        std::vector<Fxelem> factors;
        factors.push_back(pol);
        big_int r = 0;
        auto mat = formMatrixBigQ(pol);
        int n = pol.deg();
        for (int i = 0; i < n; ++i)
            mat[i][i] -= 1;
        auto base = kernelBasis<Fxelem>(mat);
        int k = base.size() + 1;//we do not have computed the first element of the base, so have to add 1 to k
        while (factors.size() < k) {
            for (int i = 0; i < factors.size(); ++i) {
                Fxelem v(base[(size_t) r]);
                for (auto &s : pol.getField().getElems()) {
                    Fxelem g = gcd(v - s, factors[i]);
                    if (g != 1 && g != factors[i]) {
                        factors[i] /= g; //We continue in the loop with the new factors[i] because it is a divisor of the old factors[i] so it is not necessary to check the previous s and r.
                        factors.push_back(g);
                        if (factors.size() == k) return factors;
                    }
                }
            }
            ++r;
        }
        return factors;
    }

    template<typename Fxelem>
    std::vector<std::pair<Fxelem, unsigned int> > factorizationBerlekamp(const Fxelem &pol) {//TODO: probar con el polinomio 1, si no funciona ponerlo como caso particular
    	std::vector<std::pair<Fxelem, unsigned int> > result;
    	if (pol.deg() <= 1){
			result.push_back(std::make_pair(pol, 1));
			return result;
		}
    	auto lc = pol.lc();
    	if (lc != 1)
    		result.push_back(Fxelem(lc), 1);
    	auto aux = squareFreeFF(pol/lc);
        for (auto &pair: aux) {
            auto aux2 = berlekamp_simple(pair.first);
            for (auto &factor: aux2) {
                result.push_back(std::make_pair(factor, pair.second));
            }
        }
        result[0].first *= lc;
        return result;
    }

    template<typename Fxelem>
    std::vector<std::pair<Fxelem, unsigned int> > factorizationCantorZassenhaus(const Fxelem &pol) {
    	std::vector<std::pair<Fxelem, unsigned int> > result;
    	if (pol.deg() <= 1){
    		result.push_back(std::make_pair(pol, 1));
    		return result;
    	}
    	auto lc = pol.lc();
    	auto aux = squareFreeFF(pol/lc);
        for (auto &pair: aux) {
            auto polAndDegree = partialFactorDD(pair.first);
            for (auto &elem: polAndDegree) {
                auto aux2 = splitFactorsDD(elem.first, elem.second);
                for (auto &factor: aux2) {
                    if (factor != 1) {
                        result.push_back(std::make_pair(factor, pair.second));
                    }
                }
            }
        }
        result[0].first *= lc;
        return result;
    }
}

#endif // __FACTORIZATION_FQ
