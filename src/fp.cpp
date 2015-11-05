#include <vector>
#include <string>
#include "fp.hpp"
#include "fpelem.hpp"
#include "exceptions.hpp"
#include "generalPurpose.hpp" // Miller Rabin

Fp::Fp(ll p): _p(p){
    // TODO create a look-up table for p < 2^16
    if(p<=0 || !millerRabin(p))
        throw EpNotPrime("Could not create F" + std::to_string(p) + ". " + std::to_string(p) + " is not prime.");
}

Fpelem Fp::get(ll n)const{return Fpelem(n, this);}

ll Fp::getSize()const{return _p;}

bool Fp::operator==(const Fp &rhs)const{return _p == rhs._p;}
bool Fp::operator!=(const Fp &rhs)const{return _p != rhs._p;}

std::vector<Fpelem> Fp::getElems(){
    std::vector<Fpelem> ret;
    for(ll i=0;i<_p;++i)
        ret.push_back(this->get(i));
    return ret;
}
