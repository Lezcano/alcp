#include "fpxelem.hpp"
#include "generalPurpose.hpp" // fastPowMod
#include "zxelem.hpp"

#include <vector>
// Fpelem included in fpxelem.hpp

// Auxiliary function for the Fpxelem ctor
Fpxelem toFpxelem(const Zxelem &e, big_int p){
    std::vector<Fpelem> v;
    auto f = Fp(p);
    for(unsigned int i=0;i<e._v.size();++i)
        v.push_back(f.get(e._v[i]));
    return v;
}

Fpxelem::Fpxelem() = default;
Fpxelem::Fpxelem(const Fpelem & e)              : PolynomialRing<Fpxelem, Fpelem>(e){}
Fpxelem::Fpxelem(const std::vector<Fpelem> & v) : PolynomialRing<Fpxelem, Fpelem>(v){}
Fpxelem::Fpxelem(const Zxelem & e, big_int p)   : PolynomialRing<Fpxelem, Fpelem>(toFpxelem(e, p)){}

const Fpxelem::F Fpxelem::getField()const{
    return this->lc().getField();
}

big_int Fpxelem::getSize()const{
    return this->getField().getSize();
}

bool Fpxelem::irreducible()const{
    Fpxelem x({this->getField().get(0), this->getField().get(1)});
    Fpxelem xpk = x; // x^(p^k)

    for(unsigned int i=0;i<this->deg()/2;++i){
        xpk = fastPowMod(xpk, this->getSize(), *this);
        if(gcd(*this, xpk-x).deg()!=0)
            return false;
    }
    return true;
}

Fpxelem getZero(const Fpxelem &e){return Fpxelem(e.getField().get(0));}
Fpxelem getOne(const Fpxelem &e){return Fpxelem(e.getField().get(1));}

const Fpelem unit(const Fpxelem &e){
	return e.lc(); }

bool compatible(const Fpxelem &lhs, const Fpxelem &rhs){
    return lhs.getField()==rhs.getField();
}

bool operator==(const Fpxelem &lhs, big_int rhs){
    return lhs.deg()==0 && lhs.lc()==lhs.getField().get(rhs);
}

bool operator==(big_int lhs, const Fpxelem &rhs){
    return rhs == lhs;
}

bool operator!=(const Fpxelem &lhs, big_int rhs){
    return !(lhs == rhs);
}

bool operator!=(big_int lhs, const Fpxelem &rhs ){
    return !(rhs == lhs);
}

