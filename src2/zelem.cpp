#include <string>
#include <vector>

#include "zelem.hpp"
#include "zxelem.hpp"
// types.hpp included en zxelem.hpp

////////////////////////////// Zelem //////////////////////////////
bool compatible(bint lhs, bint rhs){
    return true;
}
const bint unit(bint e){
    return e >= 0 ? 1 : -1;
}
const bint normalForm(bint e){ return e/unit(e);}

std::string to_string(bint e){return std::to_string(e);}
bint getZero(bint e){return 0;}
bint getOne(bint e){return 1;}

////////////////////////////// Zelem //////////////////////////////

////////////////////////////// Zxelem //////////////////////////////

Zxelem::Zxelem(const bint & e) : PolinomialRing<Zxelem, bint>(e){}
Zxelem::Zxelem(const std::vector<bint> & v) : PolinomialRing<Zxelem, bint>(v){}


const bint unit(const Zxelem &e){
    return unit(e.lc());
}

////////////////////////////// Zxelem //////////////////////////////

