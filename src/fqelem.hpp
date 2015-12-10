#ifndef __FQELEM_HPP
#define __FQELEM_HPP

#include "felembase.hpp"
#include "fq.hpp"
#include "fpxelem.hpp"
#include "types.hpp"

template<class Integer>
class Fq;

template<class Integer = big_int>
class Fqelem : public FelemBase<Fqelem<Integer>, Fq<Integer>, Fpxelem<Integer>, Integer>{
    private:
        typedef FelemBase<Fqelem, Fq<Integer>, Fpxelem<Integer>, Integer> FBase;

    public:
        using FBase::FelemBase;
        Fqelem() = default;

        Fqelem & operator=(const Fqelem &rhs){
            if(&rhs != this){
                FBase::operator=(rhs);
            }
            return *this;
        }

        Fqelem & operator=(Integer rhs){
            FBase::operator=(rhs);
            return *this;
        }

    private:
        friend class Fq<Integer>;
        Fqelem(const Fpxelem<Integer> n, Fq<Integer> f) : FBase(n, f){}
};

using Fqelem_b = Fqelem<big_int>;

#endif // __FQELEM_HPP
