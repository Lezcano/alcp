// Implementation of a GF(p) field
#ifndef __FPXELEM2_HPP
#define __FPXELEM2_HPP

#include "fp.hpp"
#include "fpelem.hpp"
#include "polRing.hpp"


class Fpxelem : public PolynomialRing<Fpxelem, Fpelem>{
    public:
        // Base field
        using F = Fp;
        using Felem = Fpelem;

        Fpxelem(const Fpelem & e);
        Fpxelem(const std::vector<Fpelem> & v);

        bool irreducible()const;
        const F getField()const;
        big_int getSize()const;
};

Fpxelem getZero(const Fpxelem &e);
Fpxelem getOne(const Fpxelem &e);
const Fpelem unit(const Fpxelem &e);
bool compatible(const Fpxelem &lhs, const Fpxelem &rhs);
bool operator==(const Fpxelem &lhs, ll rhs);
bool operator==(ll lhs, const Fpxelem &rhs);
bool operator!=(const Fpxelem &lhs, ll rhs);
bool operator!=(ll lhs, const Fpxelem &rhs );

#endif
