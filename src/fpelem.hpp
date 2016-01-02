#ifndef __FPELEM_HPP
#define __FPELEM_HPP

#include <string>
#include <vector>

#include "types.hpp"
#include "exceptions.hpp"
#include "zelem.hpp"            // to_string
#include "generalPurpose.hpp"   // millerRabin
#include "quotientRing.hpp"

namespace alcp {
    template<class Integer>
    class Fpelem;

    template<class Integer = big_int>
    class Fp {
    public:
        Fp(Integer p) : _p(p) {
            using std::to_string;
            // TODO create a look-up table for p < 2^16
            if (p <= 0 || !millerRabin(p))
                throw EpNotPrime("Could not create F" + to_string(p) + ". " + to_string(p) + " is not prime.");
        }

        Fp(const Fp<Integer> &f) = default;

        Fpelem<Integer> get(big_int n) const {
            n %= _p;
            if (n < 0)
                n += _p;
            return Fpelem<Integer>(n, *this);
        }

        Integer mod() const { return _p; }

        Integer getSize() const { return _p; }

        Integer getP() const { return _p; }

        std::size_t getM() const { return 1; }

        std::vector<Fpelem<Integer>> getElems() const {
            std::vector<Fpelem<Integer>> ret(_p);
            for (Integer i = 0; i < _p; ++i)
                ret[i] = this->get(i);
            return ret;
        }

        bool operator==(const Fp &rhs) const { return _p == rhs._p; }

        bool operator!=(const Fp &rhs) const { return _p != rhs._p; }

        friend std::string to_string(const Fp<Integer> &f) {
            using std::to_string;
            return "F" + to_string(f._p);
        }

    private:
        Integer _p;
    };


    template<class Integer>
    class Fpelem : public QuotientRing<Fpelem, Fp, Integer, Integer> {
    private:
        using FBase = QuotientRing<Fpelem, Fp, Integer, Integer>;

    public:
        using FBase::QuotientRing;
        using FBase::operator=;

        Fpelem() = default;

    private:
        friend class Fp<Integer>;

        Fpelem(const Integer n, Fp<Integer> f) : FBase(n, f) { }
    };

    using Fpelem_b = Fpelem<big_int>;
    using Fp_b = Fp<big_int>;
}

#endif // __FPELEM_HPP
