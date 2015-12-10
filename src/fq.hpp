#ifndef __FIELDQ_HPP
#define __FIELDQ_HPP

#include "types.hpp"
#include "fp.hpp"
#include "fpxelem.hpp"

#include <vector>

class Fqelem;

class Fq{
    public:
        Fq(big_int p, int n);

        Fqelem get(big_int n)const;

        Fqelem get(Fpxelem f)const;

        Fpxelem mod() const;
        big_int getSize()const;
        big_int getP()const;
        big_int getM()const;


        std::vector<Fqelem> getElems()const;

        bool operator==(const Fq &rhs)const;
        bool operator!=(const Fq &rhs)const;

        friend std::string to_string(const Fq &e);
    private:

        big_int _p;
        int _n;
        Fp _base;
        Fpxelem _mod;
};

#endif // __FIELDQ_HPP
