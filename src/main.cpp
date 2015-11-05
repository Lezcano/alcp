#include<iostream>
#include<vector>
//#include"gtest/gtest.h"
#include"fp.hpp"
#include"fpelem.hpp"
#include"fpxelem.hpp"
#include"exceptions.hpp"
#include "berlekamp.hpp"
using namespace std;

int main (){
    /*
    Fp field(5);
    Fpelem f = field.get(3);
    Fpelem g = field.get(4);
    cout << f << endl;
    cout << f+g << endl;
    cout << f + field.get(2) << endl;
    cout << f + 3 << endl;

    Fp field2(17);
    vector<Fpelem> v, w;
    for(int i=0;i<5;++i)
        v.push_back(field2.get(i+2));
    for(int i=0;i<6;++i){
        w.push_back(field2.get(3*i+3));
    }
    Fpxelem a(v), b(w);
    cout << "a\n" << a << endl;
    cout << "b\n" << b << endl;

   try{
        //cout << Fpelem(3,6)<<endl; //Throws except
        cout << std::boolalpha << (Fp(5).get(0) == Fp(5).get(5)) << endl;
        cout << "Suma a+b\n" << a+b << endl;
        cout << "suma b+a\n" << b+a << endl;
        cout << "Resta b-a\n" <<b-a << endl;
        cout << "Resta a-b\n" << a-b << endl;
        cout << "Multip a*b\n" << a*b << endl;
        std::pair<Fpxelem,Fpxelem> div1 = a.div2(b);
        std::cout << "Div & Mod a/b" << std::endl;
        std::cout << div1.second << std::endl;
        std::cout << div1.first << std::endl;

        std::pair<Fpxelem,Fpxelem> div2 = b.div2(a);
        std::cout << "Div & Mod b/a" << std::endl;
        std::cout << div2.first << std::endl << div2.second << std::endl;
        // Until here everything is ok
        */
    vector<Fpelem> v, w;
    Fp field(11);
    try{
		int c[7]= {1, -3, -1, -3, 1, -3, 1};
        for(int i=0;i<7;++i)
            v.push_back(field.get(c[i]));
		Fpxelem a(v);
		cout << a << endl;

		std::vector< Fpxelem > qw (berlekamp_simple (a));
		cout << qw[0];
		cout << qw[1];
		cout << qw[2];
		/*
        for(int i=0;i<5;++i)
            v.push_back(field.get(7*i));
        for(int i=0;i<6;++i){
            w.push_back(field.get(5*i+2));
        }
        Fpxelem a(v), b(w);
        cout << a.deg() << endl;
        cout << b.deg() << endl;
        cout << "a: " << a << endl;
        cout << "b: " << b << endl;

        cout << "Suma a+b\n" << a+b << endl;
        cout << "suma b+a\n" << b+a << endl;
        cout << "Resta b-a\n" <<b-a << endl;
        cout << "Resta a-b\n" << a-b << endl;
        cout << "Multip a*b\n" << a*b << endl;
        std::pair<Fpxelem,Fpxelem> div1 = a.div2(b);
        std::cout << "Div & Mod a/b" << std::endl;
        std::cout << div1.second << std::endl;
        std::cout << div1.first << std::endl;

        std::pair<Fpxelem,Fpxelem> div2 = b.div2(a);
        std::cout << "Div & Mod b/a" << std::endl;
        std::cout << div2.first << std::endl << div2.second << std::endl;
*/

   }catch(ExcepALCP e){
    cout << e << endl;
    }

}
/*
TEST(AddPolTest, Some){
    vector<Fpelem> v, w;
    Fp f(17);
    for(int i=0;i<5;++i)
        v.push_back(f.get(7*i));
    for(int i=0;i<6;++i){
        w.push_back(f.get(5*i+2));
    }
    Fpxelem a(v), b(w);
    EXPECT_EQ(Fpxelem({f.get(2),f.get(14),f.get(9),f.get(4),f.get(16),f.get(10)}), a+b);
}*/
