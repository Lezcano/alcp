// Implementation of a GF(p) field
#ifndef __FPXELEM2_HPP
#define __FPXELEM2_HPP

#include "types.hpp"
#include "fp.hpp"
#include "fpelem.hpp"
#include "zxelem.hpp"
#include "polRing.hpp"

class Zxelem;

class Fpxelem : public PolynomialRing<Fpxelem, Fpelem>{
    public:
        // Base field
        using F = Fp;
        using Felem = Fpelem;

        Fpxelem(const Fpelem & e);
        Fpxelem(const std::vector<Fpelem> & v);
        Fpxelem(const Zxelem & e, big_int p);

        bool irreducible()const;
        const F getField()const;
        big_int getSize()const;

        friend class Zxelem;
        friend Zxelem toZxelemSym(const Fpxelem &e);

};

Fpxelem getZero(const Fpxelem &e);
Fpxelem getOne(const Fpxelem &e);
const Fpelem unit(const Fpxelem &e);
bool compatible(const Fpxelem &lhs, const Fpxelem &rhs);
bool operator==(const Fpxelem &lhs, big_int rhs);
bool operator==(big_int lhs, const Fpxelem &rhs);
bool operator!=(const Fpxelem &lhs, big_int rhs);
bool operator!=(big_int lhs, const Fpxelem &rhs );

#endif
