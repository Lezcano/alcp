#include <vector>

#include "zxelem.hpp"
// types.hpp included en zxelem.hpp

Zxelem::Zxelem(const bint & e) : PolynomialRing<Zxelem, bint>(e){}
Zxelem::Zxelem(const std::vector<bint> & v) : PolynomialRing<Zxelem, bint>(v){}


const bint unit(const Zxelem &e){
    return unit(e.lc());
}



