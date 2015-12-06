#include <iostream>
#include <vector>
#include <stdexcept>
//#include"gtest/gtest.h"


#include "types.hpp"
#include "fp.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "fq.hpp"
#include "fqelem.hpp"
#include "zxelem.hpp"
#include "exceptions.hpp"
#include "factorizationFq.hpp"
#include "integerCRA.hpp"
#include "hensel.hpp"

using namespace std;
using namespace boost::multiprecision;

void pruebas(){
    try{
        Fp field(6);
    }
    catch(ExcepALCP& e){}


}

void symForm(){
    Fp f(17);
    vector<Fpelem> v;
    for(int i = 0; i < 9; ++i)
        v.push_back(f.get(7*i+3));
    Fpxelem aux (v);
    Zxelem zx (aux);
    Zxelem test({3,-7,0, 7, -3, 4, -6, 1, 8});
    if (test != zx)
        throw std::runtime_error("Falla el paso a forma simetrica!");
    else
        cout << "Ok. Symmetric form." << endl;
}

void pruebasHenselSqFree(){
	/*std::vector<big_int> v(17, 0);
	v[16] = 1;
	v[4] = 11;
	v[0] = 121;
	*/
	//std::vector<big_int> v = {0, 4, 22, -44, 3, 20, -2, -4, 1};
	//std::vector<big_int> v = {121, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
	//std::vector<big_int> v = {576, 0, -960, 0, 352, 0, -40, 0, 1};
	//std::vector<big_int> v = {0, 4, 6, 4, 1};
	//std::vector<big_int> v = {-24, -46, -78, -109, -115, -99, -75, -42, -11, -5, 3, 1};//(x-3) (x+4) (x^2+2) (x+1) (x^2+1) (x^4+x^3+x^2+x+1)
	std::vector<big_int> v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};//Ciclotomicos de 13 y 11
	//std::vector<big_int> v = {0, 1, 6 , 20 , 49 , 99 , 175 , 280 , 414 , 574 , 755 , 951 , 1155 , 1359 , 1554 , 1730 , 1876 , 1981 , 2036 , 2036 , 1981 , 1876 , 1730 , 1554 , 1359 , 1155 , 951 , 755 , 574 , 414 , 280 , 175 , 99 , 49 , 20 , 6 , 1};//(x^12+x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1)*(x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1)*(x^6+x^5+x^4+x^3+x^2+x+1)*(x^4+x^3+x^2+x+1)*(x^2+x+1)*(x+1)*x

	Zxelem u(0), w(0);
	Zxelem pol(v);
	//HenselLifting(pol, Fpxelem(fp), Fpxelem(ffp), u, w);
	//cout << Zxelem(r)*Zxelem(t) << endl;
	cout << pol << endl << endl;
	auto a =  factorizationHenselSquareFree(pol);
	Zxelem aux(1);
	for (auto p: a){
		aux*=p;
		cout << p << endl;
	}
	if (aux == pol){
		cout << "Bien factorizado!!" << endl;
	}
	else{
		cout << ":( Lo siento, no se ha factorizado bien." << endl;

	}

}

int main (){
	/*
	Fp field(5);
	int d[12]= {5040, -432, 10, 1};
	vector<big_int> v;
	vector<Fpelem> u1, w1;
	for(int i=0;i<=3;++i)
			v.push_back(d[i]);

	Zxelem a(v), b(field.get(0)), c(field.get(0));
	u1.push_back(field.get(0));
	u1.push_back(field.get(1));
	w1.push_back(field.get(-2));
	w1.push_back(field.get(0));
	w1.push_back(field.get(1));
	Fpxelem u(u1), w(w1);
	*/

    pruebasHenselSqFree();
//    symForm();

    /*
    // No compila
	if (HenselLifting(a, 5, u, w, b, c)){
		cout << "u: " << b << endl;
		cout << "w: " << c << endl;
	}
	else
		cout << "No such factorization" << endl;
        */
}
/*
	int c[12]= {1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1};
	Fp field(3);
	vector<Fpelem> v;
	for(int i=0;i<=11;++i)
		v.push_back(field.get(c[i]));
	vector<Fpelem> w;
	w.push_back(field.get(1));
	w.push_back(field.get(1));
	Fpxelem a(v), b(w);
	a *= b*b*b;
	auto factors = factorizationBerlekamp (a);
	for (auto &pair : factors){
			cout << "(" << pair.first << ")^" << pair.second<< endl;
	}
*/

