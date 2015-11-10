#include <string>

#include "zelem.hpp"
// Types inlcuded in zelem.hpp

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
