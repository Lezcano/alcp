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

        Fqxelem (const std::vector<std::vector<Integer>> & v, const Fq<Integer>& f) : FBase(
                [](const std::vector<std::vector<Integer>>& v1, const Fq<Integer>& f1) -> Fqxelem {
                    std::vector<Fqelem<Integer>> ret(v1.size());
                    std::transform(v1.begin(),
                                   v1.end(),
                                   ret.begin(),
                                   [&f1](const std::vector<Integer>& v2){
                                       return f1.get(Fpxelem<Integer>(v2,f1.getP()));
                                   });
                    return ret;
                }(v,f)){}

        Fqxelem (std::vector<std::vector<Integer>> && v, Fq<Integer>&& f) : FBase(
                [](std::vector<std::vector<Integer>>&& v1, Fq<Integer>&& f1) -> Fqxelem {
                    std::vector<Fqelem<Integer>> ret(v1.size());
                    std::transform(std::make_move_iterator(v1.begin()),
                                   std::make_move_iterator(v1.end()),
                                   ret.begin(),
                                   [&f1](std::vector<Integer>&& v2){
                                       return f1.get(Fpxelem<Integer>(std::move(v2),f1.getP()));
                                   });
                    return ret;
                }(v,f)){}

        const F getField() const {
            return this->lc().getField();
        }

        Integer getSize() const {
            return this->getField().getSize();
        }

        friend Fqxelem<Integer> getZero(const Fqxelem<Integer> &e) { return Fqxelem<Integer>(getZero(e.lc()); }

        friend Fqxelem<Integer> getOne(const Fqxelem<Integer> &e) { return Fqxelem<Integer>(getOne(e.lc()); }

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
