#include "types.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "fqelem.hpp"
#include "zxelem.hpp"
#include "exceptions.hpp"
#include "factorizationFq.hpp"
#include "integerCRA.hpp"
#include "hensel.hpp"
#include "modularGCD.hpp"
#include "generalPurpose.hpp"

#include <iostream>
#include <vector>
#include <stdexcept>
//#include"gtest/gtest.h"

using namespace std;
using namespace boost::multiprecision;
using namespace alcp;

void symForm(){
    Fp_b f(17);
    vector<Fpelem_b> v;
    for(int i = 0; i < 9; ++i)
        v.push_back(f.get(7*i+3));
    Fpxelem_b aux (v);
    Zxelem_b zx (aux);
    Zxelem_b test({3,-7,0, 7, -3, 4, -6, 1, 8});
    if (test != zx)
        throw runtime_error("XX. Falla el paso a forma simetrica!");
    else
        cout << "Ok. Symmetric form seems to work." << endl;
}

void testModularGCD(){
    const int n = 3;
    int i;
    Zxelem_b a[n] = {Zxelem_b({-360, -171, 145, 25, 1}),
                     Zxelem_b({-5,2,8,-3,-3,0,1,0,1}),
                     Zxelem_b(vector<big_int>({2,2}))};
    Zxelem_b b[n] = {Zxelem_b({-15,-14,-1,15,14,1}),
                     Zxelem_b({21,-9,-4,0,5,0,3}),
                     Zxelem_b({1,2,1})};
    Zxelem_b res[n] = {Zxelem_b({15,14,1}),
                       Zxelem_b({1}),
                       Zxelem_b(vector<big_int>({1,1}))};
    for(i=0;i<n;++i)
        if(res[i] != modularGCD(a[i], b[i])){
            cout << "XX. ModularGCD fails on test case " << i << "." << endl;
            break;
        }
    if(i == n)
        cout << "Ok. ModularGCD seems to work." << endl;
}

void testCRA(){
    if(integerCRA({99,97,95}, {49,-21,-30}) != -272300)
        cout << "XX. IntegerCRA fails miserably." << endl;
    else
        cout << "Ok. IntegerCRA seems to work." << endl;
}


void testGoodOldGCD(){
    int i;

    for(i=0;i<6;++i){
        if(gcd(0,i) != i || gcd(i,0) != i)
            break;
    }
    if(i != 6)
        cout << "XX. gcd(" << i << ",0) or gcd(" << i << ",0) fail" << endl;
    else
        cout << "Ok. gcd(i,0) and gcd(0,i), 0<=i<6" << endl;


    for(i=-1;i>-6;--i){
        if(gcd(0,i) != -i || gcd(i,0) != -i)
            break;
    }
    if(i != -6)
        cout << "XX. gcd(" << i << ",0) or gcd(" << i << ",0) fail" << endl;
    else
        cout << "Ok. gcd(i,0) and gcd(0,i), -6<i<0" << endl;

    // Random values
    if(gcd(42,56) != 14)
        cout << "XX. gcd(42,56) != 14" << endl;
    else
        cout << "Ok. gcd(42,56)" << endl;
    if(gcd(987,1491) != 21)
        cout << "XX. gcd(987,1491) != 21" << endl;
    else
        cout << "Ok. gcd(987,1491)" << endl;
}

void auxFactor(int a, bool verbose){
    int fact;
    fact = pollardRhoBrent(a);
    if(a % fact) {
        cout << "XX. Pollard rho factoring algorithm doesn't work with " << a << endl;
        throw std::runtime_error("Error.\n");
    }
    if(fact != 1 && fact != a){
        auxFactor(fact,verbose);
        auxFactor(a/fact, verbose);
    }
    else if(verbose)
        cout << std::max(fact, a) << " ";
}

void pollardoFactorTest(bool verbose){
    bool error = false;
    vector<int> aux{3,5,9,18,26,900};
    for(auto& a : aux){
        if(verbose)
            cout << a << " ";
        try{
            auxFactor(a, verbose);
        }catch(...){
            error = true;
        }
        if(verbose)
            cout << endl;
    }
    if(!error)
        cout << "Ok. Pollard rho factoring algorithm" << endl;
}

void pollardoLogarithm(bool verbose){
    long long log;

    // Compute log_2(5) in F_2019
    pollardRhoLogarithm(2, 5, 1019, log);
    if(verbose) {
        cout << "Log " << log << endl;
        cout << "Pow " << fastPowMod<big_int, big_int>(2, log, 1019) << endl;
    }
    if(fastPowMod<long long, long long>(2,log,1019) == 5)
        cout << "Ok. Pollard rho discrete logarithm algorithm" << endl;
    else
        cout << "XX. Pollard rho fails to compute log_2(5) in F_2019" << endl;


}
void pruebasHenselSqFree(bool verbose){
	//std::vector<big_int> v = {0, 4, 22, -44, 3, 20, -2, -4, 1};
	//std::vector<big_int> v = {121, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
	//std::vector<big_int> v = {576, 0, -960, 0, 352, 0, -40, 0, 1};
	//std::vector<big_int> v = {0, 4, 6, 4, 1};
	//std::vector<big_int> v = {-24, -46, -78, -109, -115, -99, -75, -42, -11, -5, 3, 1};//(x-3) (x+4) (x^2+2) (x+1) (x^2+1) (x^4+x^3+x^2+x+1)
	//std::vector<big_int> v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 11, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};//Ciclotomicos de 13 y 11
	std::vector<big_int> v = {0, 1, 6 , 20 , 49 , 99 , 175 , 280 , 414 , 574 , 755 , 951 , 1155 , 1359 , 1554 , 1730 , 1876 , 1981 , 2036 , 2036 , 1981 , 1876 , 1730 , 1554 , 1359 , 1155 , 951 , 755 , 574 , 414 , 280 , 175 , 99 , 49 , 20 , 6 , 1};//(x^12+x^11+x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1)*(x^10+x^9+x^8+x^7+x^6+x^5+x^4+x^3+x^2+x+1)*(x^6+x^5+x^4+x^3+x^2+x+1)*(x^4+x^3+x^2+x+1)*(x^2+x+1)*(x+1)*x

	Zxelem_b pol(v);
    if(verbose)
        cout << pol << endl << endl;
	auto a =  factorizationHenselSquareFree(pol);
	Zxelem_b aux(1);
    for (auto p : a) {
        aux *= p;
        if(verbose)
            cout << p << endl;
    }
	if (aux == pol){
		cout << "Ok. Bien factorizado!!" << endl;
	}
	else{
		cout << ":( Lo siento, no se ha factorizado bien." << endl;

	}

}

void pruebasHenselNoSeparable(){

	std::vector<big_int> v = {1, 2, 1};

	Zxelem_b u(0), w(0);
	Zxelem_b pol(v);
	cout << pol << endl << endl;
	auto a =  factorizationHensel(pol);
	Zxelem_b aux(1);
	for (auto p : a){
		for (unsigned int i =0; i< p.second; i++){
			aux*=p.first;
		}
		cout << "(" << p.first << ")^" << p.second << endl;
	}
	if (aux == pol){
		cout << "Bien factorizado!!" << endl;
	}
	else{
		cout << ":( Lo siento, no se ha factorizado bien." << endl;

	}

}

int main () {
    testCRA();
    testModularGCD();
    pollardoFactorTest(false);
    pollardoLogarithm(false);
    pruebasHenselSqFree(false);
   // pruebasHenselNoSeparable();
    testGoodOldGCD();
    symForm();
}

