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
	Fp field(11);
	vector<Fpelem> v;
	for(int i=0;i<=11;++i)
		v.push_back(field.get(c[i]));
	Fpxelem a(v);
	cout << a << endl;
	auto factors = squareFreeFF (a);
	for (auto pair : factors){
			cout << "(" << pair.first << ")^" << pair.second<< endl;
	}

}
