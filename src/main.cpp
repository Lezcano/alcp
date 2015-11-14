#include <iostream>
#include <vector>
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

int main (){
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
	cout << a << endl << u << endl << w << endl;
	if (HenselLifting(a, 5, u, w, b, c)){
		cout << "u: " << b << endl;
		cout << "w: " << c << endl;
	}
	else
		cout << "No such factorization" << endl;
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

