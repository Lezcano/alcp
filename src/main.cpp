#include<iostream>
#include"fpelem.hpp"
#include"fpxelem.hpp"

using namespace std;

int main (){
    Fpelem f(3,5);
    Fpelem g(4,5);
    cout << f << endl;
    cout << f+g << endl;
    cout << f + Fpelem(2,5) << endl;
    cout << f + 3 << endl;
    vector<Fpelem> v, w;
    for(int i=0;i<5;++i)
        v.push_back(Fpelem(i+2,17));
    for(int i=0;i<6;++i){
        w.push_back(Fpelem(3*i+3,17));
        cout << w.back() << endl;}
    Fpxelem a(v), b(w);
    cout << "a\n" << a << endl;
    cout << "b\n" << b << endl;

   try{
    //cout << Fpelem(3,6)<<endl; //Throws except
    cout << (Fpelem(0,7)==Fpelem(7,7)) << endl;
    cout << "Suma a+b\n" << a+b << endl;
    cout << "suma b+a\n" << b+a << endl;
    cout << "Resta b-a\n" <<b-a << endl;
    cout << "Resta a-b\n" << a-b << endl;
    cout << "Multip a*b\n" << a*b << endl;
    std::pair<Fpxelem,Fpxelem> div2 = b.div2(a);
    std::cout << "Div & Mod b/a" << std::endl;
    std::cout << div2.first << std::endl << div2.second << std::endl;
    std::pair<Fpxelem,Fpxelem> div1 = a.div2(b);
    std::cout << "Div & Mod a/b" << std::endl;
    std::cout << div1.first << std::endl << div1.second << std::endl;

   }catch(ExcepALCP e){
    cout << e << endl;
    }

}
