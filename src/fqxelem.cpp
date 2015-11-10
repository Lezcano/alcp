#include "fqxelem.hpp"
#include "generalPurpose.hpp" // fastPowMod

#include <vector>
// Fqelem included in fpxelem.hpp

Fqxelem::Fqxelem(const Fqelem & e) : PolinomialRing<Fqxelem, Fqelem>(e){}
Fqxelem::Fqxelem(const std::vector<Fqelem> & v) : PolinomialRing<Fqxelem, Fqelem>(v){}

const Fqxelem::F Fqxelem::getField()const{
    return this->lc().getField();
}

ll Fqxelem::getSize()const{
    return this->getField().getSize();
}

Fqxelem getZero(const Fqxelem &e){return Fqxelem(e.getField().get(0));}
Fqxelem getOne(const Fqxelem &e){return Fqxelem(e.getField().get(1));}


const Fqelem unit(const Fqxelem &e){ return e.lc(); }

bool compatible(const Fqxelem &lhs, const Fqxelem &rhs){
    return lhs.getField()==rhs.getField();
}

bool operator==(const Fqxelem &lhs, ll rhs){
    return lhs.deg()==0 && lhs.lc()==lhs.getField().get(rhs);
}

bool operator==(ll lhs, const Fqxelem &rhs){
    return rhs == lhs;
}

bool operator!=(const Fqxelem &lhs, ll rhs){
    return !(lhs == rhs);
}

bool operator!=(ll lhs, const Fqxelem &rhs ){
    return !(rhs == lhs);
}

