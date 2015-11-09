#include <vector>

#include "zxelem.hpp"
// types.hpp included en zxelem.hpp

Zxelem::Zxelem(const bint & e) : PolinomialRing<Zxelem, bint>(e){}
Zxelem::Zxelem(const std::vector<bint> & v) : PolinomialRing<Zxelem, bint>(v){}


const bint unit(const Zxelem &e){
    return e.lc() >= 0 ? 1 : -1;
}



