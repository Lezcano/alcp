// Implementation of a GF(p) field
#ifndef __FPXELEM2_HPP
#define __FPXELEM2_HPP

#include "types.hpp"
#include "fp.hpp"
#include "fpelem.hpp"
#include "zxelem.hpp"
#include "polRing.hpp"

template<class Integer> class Zxelem;

template<class Integer = big_int>
class Fpxelem : public PolynomialRing<Fpxelem<Integer>, Fpelem<Integer>>{
    private:
        using FBase = PolynomialRing<Fpxelem<Integer>, Fpelem<Integer>>;

    public:
        // Base field
        using F = Fp<Integer>;
        using Felem = Fpelem<Integer>;
        // Inherit consturctors
        using FBase::PolynomialRing;

        Fpxelem() = default;
        Fpxelem(const Zxelem<Integer> & e, Integer p) : FBase(toFpxelem(e, p)){}

        bool irreducible()const{
            Fpxelem x({this->getField().get(0), this->getField().get(1)});
            Fpxelem xpk = x; // x^(p^k)

            for(unsigned int i=0;i<this->deg()/2;++i){
                xpk = fastPowMod(xpk, this->getSize(), *this);
                if(gcd(*this, xpk-x).deg()!=0)
                    return false;
            }
            return true;
        }

        const F getField()const{
            return this->lc().getField();
        }

        Integer getSize()const{
            return this->getField().getSize();
        }

        friend class Zxelem<Integer>;
        friend Zxelem<Integer> toZxelemSym(const Fpxelem<Integer> &e);

        // non-member functions
        friend Fpxelem getZero(const Fpxelem<Integer> &e){ return Fpxelem<Integer>(e.getField().get(0)); }
        friend Fpxelem getOne(const Fpxelem<Integer> &e){ return Fpxelem<Integer>(e.getField().get(1)); }
        friend const Felem unit(const Fpxelem<Integer> &e){ return e.lc(); }
        friend bool compatible(const Fpxelem<Integer> &lhs, const Fpxelem<Integer> &rhs){
            return lhs.getField()==rhs.getField();
        }
        friend bool operator==(const Fpxelem<Integer> &lhs, Integer rhs){
            return lhs.deg()==0 && lhs.lc()==lhs.getField().get(rhs);
        }
        friend bool operator==(Integer lhs, const Fpxelem<Integer> &rhs){
            return rhs == lhs;
        }
        friend bool operator!=(const Fpxelem<Integer> &lhs, Integer rhs){
            return !(lhs == rhs);
        }
        friend bool operator!=(Integer lhs, const Fpxelem<Integer> &rhs ){
            return !(rhs == lhs);
        }

    private:

        Fpxelem toFpxelem(const Zxelem<Integer> &e, Integer p){
            std::vector<Felem> v(e.deg()+1);
            auto f = Fp<Integer>(p);
            for(std::size_t i=0;i<=e.deg();++i)
                v[i] = f.get(e[i]);
            return v;
        }
};

using Fpxelem_b = Fpxelem<big_int>;

#endif
