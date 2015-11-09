#ifndef __FIELDQ_HPP
#define __FIELDQ_HPP

#include <vector>

#include "types.hpp"

class Fqelem;

class Fq{
    public:
        Fq(ll p);

        Fqelem get(ll p, int n)const;

        ll getSize()const;

        std::vector<Fqelem> getElems()const;

        bool operator==(const Fq &rhs)const;
        bool operator!=(const Fq &rhs)const;

        friend std::string to_string(const Fq &e);
    private:
        bool increase(std::vector<Fpelem &act);

        ll _p;
        int _n;
        Fpelem _base
        Fpxelem _mod;
};

#endif // __FP_HPP
