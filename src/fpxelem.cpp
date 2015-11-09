#include <vector>

#include "fpxelem.hpp"
#include "generalPurpose.hpp" // fastPowMod
// Fpelem included in fpxelem.hpp

Fpxelem::Fpxelem(const Fpelem & e) : PolinomialRing<Fpxelem, Fpelem>(e){}
Fpxelem::Fpxelem(const std::vector<Fpelem> & v) : PolinomialRing<Fpxelem, Fpelem>(v){}

const Fpxelem::F Fpxelem::getField()const{
    return this->lc().getField();
}

ll Fpxelem::getSize()const{
    return this->getField().getSize();
}

// TODO comentar!!
// TODO probar!!
bool Fpxelem::irreducible()const{
    Fpxelem x({this->getField.get(0), this->getField.get(1)});
    Fpxelem xpk = x; // x^(p^k)

    for(int i=0;i<this->degree()/2){
        xpk = fastPowMod(xpk, this->getSize(), *this);
        g = gcd(*this, xpk-x);
        if(g.deg()!=0)
            return false;
    }
    return true;
}


const Fpelem unit(const Fpxelem &e){ return e.lc(); }

bool compatible(const Fpxelem &lhs, const Fpxelem &rhs){
    return lhs.getField()==rhs.getField();
}

bool operator==(const Fpxelem &lhs, ll rhs){
    return lhs.deg()==0 && lhs.lc()==lhs.getField().get(rhs);
}

bool operator==(ll lhs, const Fpxelem &rhs){
    return rhs == lhs;
}

bool operator!=(const Fpxelem &lhs, ll rhs){
    return !(lhs == rhs);
}

bool operator!=(ll lhs, const Fpxelem &rhs ){
    return !(rhs == lhs);
}

