#ifndef __FQELEM_HPP
#define __FQELEM_HPP

#include <type_traits>
#include <vector>

#include "quotientRing.hpp"
#include "zelem.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "generalPurpose.hpp"
#include "types.hpp"

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
            while (!Fpxelem<Integer>(v).irreducible()) {
                this->increment(v);
                // It's divisible by the polynomial x
                if (v[0] == 0)
                    this->increment(v);
            }
            _mod = Fpxelem<Integer>(v);
        }

        Fq(const Fq<Integer> &f) = default;

        Fqelem<Integer> get(Integer n) const {
            return Fqelem<Integer>(Fpxelem<Integer>(_base.get(n)), _mod);
        }

        Fqelem<Integer> get(Fpxelem<Integer> f) const {
            return Fqelem<Integer>(f, _mod);
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

        friend std::string to_string(const Fq<Integer> &e) {
            using std::to_string;
            using alcp::to_string;
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

        friend class Fqelem<Integer>;
        Fq(Integer p, std::size_t m, const Fpxelem<Integer>& mod) :
                _p(p), _m(m), _base(_mod[0].getField()), _mod(mod) { }

        Integer _p;
        std::size_t _m;
        Fp<Integer> _base;
        Fpxelem<Integer> _mod;
    };


    template<class Integer = big_int>
    class Fqelem : public QuotientRing<Fqelem, Fpxelem<Integer>, Integer> {
    private:
        using FBase = QuotientRing<Fqelem, Fpxelem<Integer>, Integer>;

    public:
        using FBase::QuotientRing;
        using FBase::operator=;
        using F = Fq<Integer>;

        F getField() const{ return F(p(), m(), this->_mod); }

        friend std::string to_string_coef(const Fqelem& e){
            using std::to_string;
            using alcp::to_string;
            std::string s(to_string(e));
            // Not-so-nifty hack for the sake of readablity
            if(s.find('x') == std::string::npos)
                return "+" + to_string(e._num.lc());
            return "+(" + s + ")" ;
        }

    private:
        Integer p() const{ return this->_mod.getSize(); }
        std::size_t m() const{ return this->_mod.deg(); }

        friend class Fq<Integer>;
    };

    using Fqelem_b = Fqelem<big_int>;
    using Fq_b = Fq<big_int>;
}
#endif // __FQELEM_HPP
