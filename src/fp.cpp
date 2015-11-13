#include "fp.hpp"
// types.hpp included in fp.hpp
#include "fpelem.hpp"
#include "exceptions.hpp"
#include "zelem.hpp"
#include "generalPurpose.hpp" // Miller Rabin

#include <vector>
#include <string>

Fp::Fp(big_int p): _p(p){
    // TODO create a look-up table for p < 2^16
    if(p<=0 || !millerRabin(p))
        throw EpNotPrime("Could not create F" + to_string(p) + ". " + to_string(p) + " is not prime.");
}

Fpelem Fp::get(big_int n)const{
    return Fpelem(n, std::unique_ptr<Fp>(new Fp(*this)));
}

big_int Fp::getSize()const{return _p;}

bool Fp::operator==(const Fp &rhs)const{return _p == rhs._p;}
bool Fp::operator!=(const Fp &rhs)const{return _p != rhs._p;}

std::vector<Fpelem> Fp::getElems()const{
    std::vector<Fpelem> ret;
    for(big_int i=0;i<_p;++i)
        ret.push_back(this->get(i));
    return ret;
}
