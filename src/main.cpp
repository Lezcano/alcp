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

}
