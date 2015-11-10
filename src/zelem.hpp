#ifndef __ZELEM_HPP
#define __ZELEM_HPP

#include <string>
#include "types.hpp"

std::string to_string(bint e);
bool compatible(bint lhs, bint rhs);
const bint unit(bint e);
const bint normalForm(bint e);
bint getZero(bint e);
bint getOne(bint e);

#endif // __ZELEM_HPP
