// Implementation of a GF(p) field
#ifndef __FPXELEM2_HPP
#define __FPXELEM2_HPP

#include "fpelem.hpp"
#include "fp.hpp"

class Fpxelem : PolinomialField<Fpxelem>{
    public:
        // Base field
        using F = Fp;
        using Felem = Fpelem;

}

bool compatible(Fpxelem, Fpxelem);

#endif
