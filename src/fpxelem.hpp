// Implementation of a GF(p) field
#ifndef __FPXELEM_HPP
#define __FPXELEM_HPP

#include <algorithm>    // std::transform
#include <vector>

#include "types.hpp"
#include "fpelem.hpp"
#include "zxelem.hpp"
#include "polRing.hpp"

namespace alcp {
    template<class Integer>
    class Zxelem;

    template<class Integer>
    class Fpxelem : public PolynomialRing<Fpxelem, Fpelem<Integer>, Integer> {
    private:
        using FBase = PolynomialRing<Fpxelem, Fpelem<Integer>, Integer>;

    public:
        // Base field
        using F = Fp<Integer>;
        using Felem = Fpelem<Integer>;

        // Inherit ctors
        using FBase::PolynomialRing;

        Fpxelem() = default;

        Fpxelem(const Zxelem<Integer> &e, Integer p) : FBase(toFpxelem(e, p)) { }

        bool irreducible() const {
            Fpxelem x(std::vector<Fpelem<Integer>>{getZero(this->lc()), getOne(this->lc())});
            Fpxelem xpk = x; // x^(p^k)

            for (unsigned int i = 0; i < this->deg() / 2; ++i) {
                xpk = fastPowMod(xpk, this->getSize(), *this);
                if (gcd(*this, xpk - x).deg() != 0)
                    return false;
            }
            return true;
        }

        const Fp<Integer> getField() const {
            return this->lc().getField();
        }

        Integer getSize() const {
            return this->getField().getSize();
        }

        // It returns the symmetric representation of f \in Fp[X]
        //  as an element of Z[X]
        friend Zxelem<Integer> toZxelemSym(const Fpxelem &e) {
            std::vector<Integer> v(e.deg() + 1);
            Integer p = e.getSize();
            // When p = 2, p2 = 1 and so every v[i] is <= p2
            // When p is odd, p/2 = (p-1)/2 and so we get the
            //  symmetric representation
            Integer p2 = p / 2;
            std::transform(e._v.begin(), e._v.end(), v.begin(),
                           [&p, &p2](const Fpelem<Integer> &a) -> Integer {
                               return static_cast<Integer>(a) <= p2 ?
                                      static_cast<Integer>(a) :
                                      static_cast<Integer>(a) - static_cast<Integer>(p);
                           });
            return v;
        }


        // non-member functions
        friend Fpxelem getZero(const Fpxelem<Integer> &e) { return Fpxelem<Integer>(getZero(e.lc())); }

        friend Fpxelem getOne(const Fpxelem<Integer> &e) { return Fpxelem<Integer>(getOne(e.lc())); }

        friend Fpxelem unit(const Fpxelem<Integer> &e) { return e.lc(); }

        friend bool compatible(const Fpxelem<Integer> &lhs, const Fpxelem<Integer> &rhs) {
            return lhs.getField() == rhs.getField();
        }

        friend bool operator==(const Fpxelem<Integer> &lhs, Integer rhs) {
            return lhs.deg() == 0 && lhs.lc() == lhs.getField().get(rhs);
        }

        friend bool operator==(Integer lhs, const Fpxelem<Integer> &rhs) {
            return rhs == lhs;
        }

        friend bool operator!=(const Fpxelem<Integer> &lhs, Integer rhs) {
            return !(lhs == rhs);
        }

        friend bool operator!=(Integer lhs, const Fpxelem<Integer> &rhs) {
            return !(rhs == lhs);
        }
    };

    using Fpxelem_b = Fpxelem<big_int>;

    //template class Fpxelem<big_int>;
}

#endif
