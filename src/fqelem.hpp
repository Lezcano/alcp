#ifndef __FQELEM_HPP
#define __FQELEM_HPP

#include "quotientRing.hpp"
#include "fpxelem.hpp"
#include "types.hpp"
#include <type_traits>

namespace alcp {
    template<class Integer>
    class Fqelem;

    template<class Integer = big_int>
    class Fq {
    public:
        Fq(Integer p, std::size_t m) : _p(p), _m(m), _base(p) {
            std::vector<Fpelem<Integer>> v(_m + 1, _base.get(0));

            v.back() = _base.get(1);
            v[0] = _base.get(1);
            // Hopefully this ends in a reasonable amount of time
            //  Maybe set a timer?
            // The probability of getting an irreducible polynomial
            //  in each iteration is 1/n
            while (!Fpxelem<Integer>(v).irreducible()) {
                this->increment(v);
                if (v[0] == 0) // It's divisible by the polynomial p(x) = x
                    this->increment(v);
            }
            _mod = Fpxelem<Integer>(v);
        }

        Fq(const Fq<Integer> &f) = default;

        Fqelem<Integer> get(Integer n) const {
            return Fqelem<Integer>(Fpxelem<Integer>(_base.get(n)), Fq(*this));
        }

        Fqelem<Integer> get(Fpxelem<Integer> f) const {
            return Fqelem<Integer>(f, Fq(*this));
        }

        Fpxelem<Integer> mod() const { return _mod; }

        Integer getSize() const { return fastPow(_p, _m); }

        Integer getP() const { return _p; }

        std::size_t getM() const { return _m; }


        std::vector<Fqelem<Integer>> getElems() const {
            std::vector<Fqelem<Integer>> ret(this->getSize());
            std::vector<Fpelem<Integer>> act(_m, _base.get(0));
            Integer i = 0;

            do {
                ret[i++] = this->get(Fpxelem<Integer>(act));
            } while (this->increment(act));

            return ret;
        }

        bool operator==(const Fq &rhs) const { return _p == rhs._p && _m == rhs._m; }

        bool operator!=(const Fq &rhs) const { return _p != rhs._p && _m != rhs._m; }

        const F getField() const{ return F(p, m); }

        friend std::string to_string(const Fq<Integer> &e) {
            using std::to_string;
            return "F" + to_string(e.getP()) + "^" + to_string(e.getM());
        }

    private:
        // Auxiliary functions
        bool increment(std::vector<Fpelem<Integer>> &act) const {
            for (auto e = act.begin(); e != act.end(); ++e) {
                *e += 1;
                if (*e != 0)
                    return true;
            }
            return false; // We are done
        }

        Integer _p;
        std::size_t _m;
        Fp<Integer> _base;
        Fpxelem<Integer> _mod;
    };


    template<class Integer = big_int>
    class Fqelem : public QuotientRing<Fqelem, Fq, Fpxelem<Integer>, Integer> {
    private:
        using FBase = QuotientRing<Fqelem, Fq, Fpxelem<Integer>, Integer>;

    public:
        using FBase::QuotientRing;
        using FBase::operator=;

        Fqelem() = default;

    private:
        friend class Fq<Integer>;

        Fqelem(const Fpxelem<Integer> n, Fq<Integer> f) : FBase(n, f) { }
    };

    using Fqelem_b = Fqelem<big_int>;
    using Fq_b = Fq<big_int>;
}
#endif // __FQELEM_HPP
