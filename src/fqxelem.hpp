// Implementation of a GF(p) field
#ifndef __FQXELEM_HPP
#define __FQXELEM_HPP

#include <vector>

#include "fqelem.hpp"
#include "polRing.hpp"

namespace alcp {
    template<class Integer>
    class Fqxelem : public PolynomialRing<Fqxelem<Integer>, Fqelem<Integer>> {
    private:
        using FBase = PolynomialRing<Fqxelem<Integer>, Fqelem<Integer>>;
    public:
        // Base field
        using F = Fq<Integer>;
        using Felem = Fqelem<Integer>;

        // Inherit ctor
        using FBase::PolynomialRing;

        Fqxelem() = default;

        const F getField() const {
            return this->lc().getField();
        }

        Integer getSize() const {
            return this->getField().getSize();
        }

        friend Fqxelem<Integer> getZero(const Fqxelem<Integer> &e) { return Fqxelem<Integer>(e.getField().get(0)); }

        friend Fqxelem<Integer> getOne(const Fqxelem<Integer> &e) { return Fqxelem<Integer>(e.getField().get(0)); }

        friend Fqxelem<Integer> unit(const Fqxelem<Integer> &e) { return e.lc(); }

        friend bool compatible(const Fqxelem<Integer> &lhs, const Fqxelem<Integer> &rhs) {
            return lhs.getField() == rhs.getField();
        }

        friend bool operator==(const Fqxelem<Integer> &lhs, Integer rhs) {
            return lhs.deg() == 0 && lhs.lc() == lhs.getField().get(rhs);
        }

        friend bool operator==(Integer lhs, const Fqxelem<Integer> &rhs) { return rhs == lhs; }

        friend bool operator!=(const Fqxelem<Integer> &lhs, Integer rhs) { return !(lhs == rhs); }

        friend bool operator!=(Integer lhs, const Fqxelem<Integer> &rhs) { return !(rhs == lhs); }
    };

    using Fqxelem_b = Fqxelem<big_int>;
}

#endif
