#ifndef __FP_HPP
#define __FP_HPP

#include <vector>
#include "types.hpp"
#include "fpelem.hpp"

class Fpelem;

class Fp{
    public:
        Fp(ll p);

        Fpelem get(ll n)const;

        ll getSize()const;

        std::vector<Fpelem> getElems()const;

        bool operator==(const Fp &rhs)const;
        bool operator!=(const Fp &rhs)const;


    private:
        ll _p;
};

#endif // __FP_HPP
