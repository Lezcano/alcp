#ifndef __FQELEM_HPP
#define __FQELEM_HPP

#include <type_traits>
#include <vector>
#include <stdexcept>

#include "quotientRing.hpp"
#include "zelem.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "generalPurpose.hpp"
#include "types.hpp"
#include "exceptions.hpp"

namespace alcp {
    template<class Integer>
    class Fqelem;

    template<class Integer = big_int>
    class Fq {
    public:
        Fq(Integer p, std::size_t m){
            Fp_b f (p);
            std::vector<Fpelem<Integer>> v(m + 1, f.get(0));

            v.back() = f.get(1);
            v[0] = f.get(1);
            while (!Fpxelem<Integer>(v).irreducible()) {
                this->increment(v);
                // It's divisible by the polynomial x
                if (v[0] == 0)
                    this->increment(v);
            }
            _mod = Fpxelem<Integer>(v);
        }

        Fq(Fpxelem<Integer> mod) :_mod(mod){
            if(!mod.irreducible())
                throw EFPXNotIrreducible("The polinomial provided to Fq was not irreducible.");
        }


        Fq(const Fq<Integer> &f) = default;

        Fqelem<Integer> get(Integer n) const {
            return Fqelem<Integer>(Fpxelem<Integer>(this->_mod.getField().get(n)), _mod);
        }

        Fqelem<Integer> get(const Fpxelem<Integer>& f) const {
            return Fqelem<Integer>(f%_mod, _mod);
        }

        Fqelem<Integer> get(Fpxelem<Integer>&& f) const {
            return Fqelem<Integer>(std::move(f % _mod), _mod);
        }

        Fpxelem<Integer> mod() const { return _mod; }

        Integer getSize() const { return fastPow(this->getP(), this->getM()); }

        Integer getP() const { return _mod.lc().getSize(); }

        std::size_t getM() const { return _mod.deg(); }


        std::vector<Fqelem<Integer>> getElems() const {
            std::vector<Fqelem<Integer>> ret(static_cast<std::size_t>(this->getSize()));
            std::vector<Fpelem<Integer>> act(this->getM(), getZero(_mod.lc()));
            std::size_t i = 0;

            do {
                ret[i++] = this->get(Fpxelem<Integer>(act));
            } while (this->increment(act));

            return ret;
        }

        bool operator==(const Fq &rhs) const { return this->_mod == rhs._mod; }

        bool operator!=(const Fq &rhs) const { return this->_mod != rhs._mod; }

        friend std::string to_string(const Fq<Integer> &e) {
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

        // The same not-so-cool trick
        Fq(const Fpxelem<Integer>& mod, bool) : _mod(mod) { }

        Fpxelem<Integer> _mod;
    };


    template<class Integer = big_int>
    class Fqelem : public QuotientRing<Fqelem, Fpxelem<Integer>, Integer> {
    private:
        using FBase = QuotientRing<Fqelem, Fpxelem<Integer>, Integer>;

    public:
        using F = Fq<Integer>;
        using FBase::QuotientRing;
        using FBase::operator=;

        F getField() const{ return F(this->_mod, true); }

        friend std::string to_string(const Fqelem &e) {
            return to_string(e._num, 't');
        }

        friend std::string to_string_coef(const Fqelem& e){
            if(e._num.deg() == 0)
                return "+" + to_string(e._num.lc());
            else if(e._num.nonZeroCoefs() == 1)
                return "+" + to_string(e._num, 't');
            return "+(" + to_string(e._num, 't') + ")" ;
        }

    private:
        friend class Fq<Integer>;
    };

    using Fqelem_b = Fqelem<big_int>;
    using Fq_b = Fq<big_int>;
}
#endif // __FQELEM_HPP
