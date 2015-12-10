#ifndef __FPELEM_HPP
#define __FPELEM_HPP

#include "felembase.hpp"
#include "fp.hpp"
#include "zelem.hpp"
#include "types.hpp"

template<class Integer> class Fp;

template<class Integer>
class Fpelem : public FelemBase<Fpelem<Integer>, Fp<Integer>, Integer, Integer>{
    private:
        typedef FelemBase<Fpelem, Fp<Integer>, Integer, Integer> FBase;

    public:
        using FBase::FelemBase;
        Fpelem() = default;

        // I do not know why do I have to redeclare these...
        Fpelem & operator=(const Fpelem &rhs){
            if(&rhs != this){
                FBase::operator=(rhs);
            }
            return *this;
        }

        Fpelem & operator=(Integer rhs){
            FBase::operator=(rhs);
            return *this;
        }



    private:
        friend class Fp<Integer>;
        Fpelem(const Integer n, Fp<Integer> f) : FBase(n, f){}
};

using Fpelem_b = Fpelem<big_int>;

#endif // __FPELEM_HPP
