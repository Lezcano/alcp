#ifndef __ZELEM_HPP
#define __ZELEM_HPP

#include <string>
#include "types.hpp"

std::string to_string(big_int e);
bool compatible(big_int lhs, big_int rhs);
const big_int unit(big_int e);
const big_int normalForm(big_int e);
big_int getZero(big_int e);
big_int getOne(big_int e);

#endif // __ZELEM_HPP
