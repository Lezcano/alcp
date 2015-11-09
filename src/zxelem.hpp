// Implementation of a GF(p) field
#ifndef __ZXELEM_HPP
#define __ZXELEM_HPP

#include "types.hpp"
#include "zelem.hpp"    // Auxliary functions
#include "polRing.hpp"


class Zxelem : public PolinomialRing<Zxelem, bint>{
    public:
        Zxelem(const bint & e);
        Zxelem(const std::vector<bint> & v);
};

const bint unit(const Zxelem &e);

#endif
