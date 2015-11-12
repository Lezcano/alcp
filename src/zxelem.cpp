#include <vector>

#include "zxelem.hpp"
// types.hpp included en zxelem.hpp

Zxelem::Zxelem(const big_int & e) : PolynomialRing<Zxelem, big_int>(e){}
Zxelem::Zxelem(const std::vector<big_int> & v) : PolynomialRing<Zxelem, big_int>(v){}


const big_int unit(const Zxelem &e){
    return unit(e.lc());
}



