#ifndef __FPELEM_HPP
#define __FPELEM_HPP

#include <string>
#include <vector>
#include <cctype>

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
    static_assert(is_integral<Integer>::value, "Type is not a supported integer.");

    public:
        Fp(Integer p) : _p(p) {
            if (_p <= 0 || !millerRabin(_p))
                throw EpNotPrime("Could not create F" + to_string(_p) + ". " + to_string(_p) + " is not prime.");
        }

        Fp(const Fp<Integer> &f) = default;

        Fpelem<Integer> get(Integer n) const {
            n %= _p;
            if (n < 0)
                n += _p;
            return Fpelem<Integer>(n, _p);
        }


        Integer mod() const { return _p; }

        Integer getSize() const { return _p; }

        Integer getP() const { return _p; }

        std::size_t getM() const { return 1; }

        std::vector<Fpelem<Integer>> getElems() const {
            std::vector<Fpelem<Integer>> ret(static_cast<std::size_t>(_p));
            for (std::size_t i = 0; i < _p; ++i)
                ret[i] = this->get(i);
            return ret;
        }



        bool operator==(const Fp &rhs) const { return _p == rhs._p; }

        bool operator!=(const Fp &rhs) const { return _p != rhs._p; }

        friend std::string to_string(const Fp<Integer> &f) {
            return "F" + to_string(f._p);
        }

        friend std::ostream &operator<<(std::ostream &os, Fp<Integer> e) {
            os << to_string(e);
            return os;
        }

    private:
        Integer _p;

        friend class Fpelem<Integer>;

        // Small hack to let Fpelem create an Fp fast without memoization
        Fp(Integer p, bool) : _p(p) { }

    };


    template<class Integer>
    class Fpelem : public QuotientRing<Fpelem, Integer, Integer> {
    private:
        // ::alcp::Fpelem still necessary for clang 3.9
        // http://stackoverflow.com/questions/17687459/clang-not-accepting-use-of-template-template-parameter-when-using-crtp
        using FBase = QuotientRing<::alcp::Fpelem, Integer, Integer>;

    public:
        using F = Fp<Integer>;
        using FBase::FBase;
        using FBase::operator=;

        Fpelem() = default;
        Fpelem(const Fpelem & e)  = default;
        Fpelem(Fpelem && e)  = default;
        Fpelem & operator=(const Fpelem & e)  = default;
        Fpelem & operator=(Fpelem && e)  = default;

        F getField() const{ return F(this->_mod, true); }

        friend std::string to_string(const Fpelem &e) {
            return to_string(e._num);
        }

        friend std::string to_string_coef(const Fpelem& e){
            return "+" + to_string(e);
        }

    private:
        template <class>
            friend class Fp;
    };

    using Fpelem_b = Fpelem<big_int>;
    using Fp_b = Fp<big_int>;

    //template class Fpelem<big_int>;
    //template class Fp<big_int>;
}

#endif // __FPELEM_HPP
