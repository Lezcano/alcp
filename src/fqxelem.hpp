// Implementation of a GF(p) field
#ifndef __FQXELEM_HPP
#define __FQXELEM_HPP

#include "fq.hpp"
#include "fqelem.hpp"
#include "polRing.hpp"


class Fqxelem : public PolynomialRing<Fqxelem, Fqelem>{
    public:
        // Base field
        using F = Fq;
        using Felem = Fqelem;

        Fqxelem(const Fqelem & e);
        Fqxelem(const std::vector<Fqelem> & v);

        const F getField()const;
        ll getSize()const;
};

Fqxelem getZero(const Fqxelem &e);
Fqxelem getOne(const Fqxelem &e);
const Fqelem unit(const Fqxelem &e);
bool compatible(const Fqxelem &lhs, const Fqxelem &rhs);
bool operator==(const Fqxelem &lhs, ll rhs);
bool operator==(ll lhs, const Fqxelem &rhs);
bool operator!=(const Fqxelem &lhs, ll rhs);
bool operator!=(ll lhs, const Fqxelem &rhs );

#endif
