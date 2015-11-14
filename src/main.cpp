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
using namespace std;
using namespace boost::multiprecision;

int main (){

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

}
