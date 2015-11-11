#include "fpxelem.hpp"
#include "generalPurpose.hpp" // fastPowMod

#include <vector>
// Fpelem included in fpxelem.hpp

Fpxelem::Fpxelem(const Fpelem & e) : PolynomialRing<Fpxelem, Fpelem>(e){}
Fpxelem::Fpxelem(const std::vector<Fpelem> & v) : PolynomialRing<Fpxelem, Fpelem>(v){}

const Fpxelem::F Fpxelem::getField()const{
    return this->lc().getField();
}

ll Fpxelem::getSize()const{
    return this->getField().getSize();
}

bool Fpxelem::irreducible()const{
    Fpxelem x({this->getField().get(0), this->getField().get(1)});
    Fpxelem xpk = x; // x^(p^k)

    for(int i=0;i<this->deg()/2;++i){
        xpk = fastPowMod<Fpxelem>(xpk, this->getSize(), *this);
        if(gcd(*this, xpk-x).deg()!=0)
            return false;
    }
    return true;
}

Fpxelem getZero(const Fpxelem &e){return Fpxelem(e.getField().get(0));}
Fpxelem getOne(const Fpxelem &e){return Fpxelem(e.getField().get(1));}

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

