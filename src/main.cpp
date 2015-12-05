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
	std::vector<big_int> v = {34, 4, 22, -44, 3, 20, -2, -4, 1};
	std::vector<big_int> r = { 2, 4, 6, 4, 1};
	std::vector<big_int> t = { 17, -32, 24, -8, 1};
	std::vector<big_int> o = { 2, 4, 1, 4, 1};
	std::vector<big_int> p = { 2, 3, 4, 2, 1};
	std::vector<Fpelem> fp, ffp;
	for (int i =0; i<5; i++){
		fp.push_back(Fp(5).get(o[i]));
		ffp.push_back(Fp(5).get(p[i]));
	}
	Zxelem u(0), w(0);
	Zxelem pol(v);
	HenselLifting(pol, Fpxelem(fp), Fpxelem(ffp), u, w);
	/*Zxelem pol(v);
	cout << Zxelem(r)*Zxelem(t) << endl;
	auto a =  factorizationHenselSquareFree(pol);

	for (auto p: a){
		cout << p << endl;
	}
	*/
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

//    pruebasHenselSqFree();
    symForm();

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

