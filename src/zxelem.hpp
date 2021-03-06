#ifndef __ZXELEM_HPP
#define __ZXELEM_HPP

#include <vector>

#include "types.hpp"
#include "zelem.hpp"    // Auxliary functions
#include "fpxelem.hpp"
#include "polRing.hpp"

namespace alcp {
    template<class Integer>
    class Fpxelem;

    template<class Integer>
    class Fpelem;

    template<class Integer>
    class Fp;

    template<class Integer>
    class Zxelem : public PolynomialRing<Zxelem, Integer, Integer> {
    private:
        // ::alcp::Zxelem still necessary for clang 3.9
        // http://stackoverflow.com/questions/17687459/clang-not-accepting-use-of-template-template-parameter-when-using-crtp
        using RBase = PolynomialRing<::alcp::Zxelem, Integer, Integer>;

    public:
        // Inherit ctors
        using RBase::RBase;

        Zxelem() = default;

        /* It returns a polynomial in Zx using the symmetric
         * representation of Fp, this is
         * -(p-1)/2, -(p-1)/2+1, ..-1, 0, 1, .. (p-1)/2
         * */
        Zxelem(const Fpxelem<Integer> &e) : RBase(toZxelemSym(e)) { }

        friend Zxelem<Integer> getZero(const Zxelem<Integer> &e) { return Zxelem<Integer>(0); }

        friend Zxelem<Integer> getOne(const Zxelem<Integer> &e) { return Zxelem<Integer>(1); }

        friend Zxelem<Integer> unit(const Zxelem<Integer> &e) { return unit(e.lc()); }

        friend bool compatible(const Zxelem<Integer> &, const Zxelem<Integer> &) {
            return true;
        }

        friend Integer normInf(const Zxelem<Integer> &e) {
            Integer ret = e[0];
            for (std::size_t i = 1; i <= e.deg(); ++i)
                if (e[i] > ret)
                    ret = e[i];
            return ret;
        }

        friend Integer content(const Zxelem<Integer> &e) {
            Integer gcdE = e[0];
            for (size_t i = 1; i <= e.deg(); ++i)
                gcdE = gcd(gcdE, e[i]);
            return gcdE;
        }
    };

    using Zxelem_b = Zxelem<big_int>;

    //template class Zxelem<big_int>;
}
#endif
