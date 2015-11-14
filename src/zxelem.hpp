#ifndef __ZXELEM_HPP
#define __ZXELEM_HPP

#include "types.hpp"
#include "zelem.hpp"    // Auxliary functions
#include "polRing.hpp"

class Fpxelem;


class Zxelem : public PolynomialRing<Zxelem, big_int>{
    public:
        Zxelem(const big_int & e);
        Zxelem(const std::vector<big_int> & v);
        Zxelem(const Fpxelem & e);

};

Zxelem getZero(const Zxelem &e);
Zxelem getOne(const Zxelem &e);
const big_int unit(const Zxelem &e);
bool compatible(const Zxelem &lhs, const Zxelem &rhs);
big_int normInf(const Zxelem &e);

#endif
