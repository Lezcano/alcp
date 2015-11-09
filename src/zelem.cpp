#include <string>

#include "zelem.hpp"
// Types inlcuded in zelem.hpp

bool compatible(bint lhs, bint rhs){
    return true;
}
std::string to_string(bint e){return std::to_string(e);}
bint getZero(bint e){return 0;}
bint getOne(bint e){return 1;}
