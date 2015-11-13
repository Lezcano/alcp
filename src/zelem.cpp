#include <string>

#include "zelem.hpp"
// Types inlcuded in zelem.hpp

bool compatible(big_int lhs, big_int rhs){
    return true;
}
const big_int unit(big_int e){
    return e >= 0 ? (big_int)1 : (big_int)-1;
}
const big_int normalForm(big_int e){ return e/unit(e);}

std::string to_string(big_int e){return std::to_string(e);}

big_int getZero(big_int e){return (big_int)0;}
big_int getOne(big_int e){return (big_int)1;}
