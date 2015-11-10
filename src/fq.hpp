#ifndef __FIELDQ_HPP
#define __FIELDQ_HPP

#include "types.hpp"
#include "fp.hpp"
#include "fpxelem.hpp"

#include <vector>

class Fqelem;

class Fq{
    public:
        Fq(ll p, int n);

        Fqelem get(ll n)const;

        Fqelem get(Fpxelem f)const;

        ll getSize()const;

        std::vector<Fqelem> getElems()const;

        bool operator==(const Fq &rhs)const;
        bool operator!=(const Fq &rhs)const;

        friend std::string to_string(const Fq &e);
    private:

        ll _p;
        int _n;
        Fp _base;
        Fpxelem _mod;
};

#endif // __FP_HPP
