#ifndef __ZXELEM_HPP
#define __ZXELEM_HPP

#include "types.hpp"
#include "zelem.hpp"    // Auxliary functions
#include "polRing.hpp"


class Zxelem : public PolynomialRing<Zxelem, big_int>{
    public:
        Zxelem(const big_int & e);
        Zxelem(const std::vector<big_int> & v);
};

const big_int unit(const Zxelem &e);

#endif
