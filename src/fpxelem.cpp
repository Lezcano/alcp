#include <vector>
#include "fpelem.hpp"
#include "fpxelem.hpp"

Fpxelem::Fpxelem(const Fpelem & e) : PolinomialRing<Fpxelem, Fpelem>(e){}
Fpxelem::Fpxelem(const std::vector<Fpelem> & v) : PolinomialRing<Fpxelem, Fpelem>(v){}

const Fpxelem::F Fpxelem::getField()const{
    return this->lc().getField();
}

ll Fpxelem::getSize()const{
    return this->getField().getSize();
}


const Fpelem unit(const Fpxelem &e){ return e.lc(); }
const Fpxelem normalForm(const Fpxelem &e){ return e/unit(e); }

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

