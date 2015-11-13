#ifndef __FIELDP_HPP
#define __FIELDP_HPP

#include "types.hpp"

#include <vector>

class Fpelem;

class Fp{
    public:
        Fp(big_int p);

        Fpelem get(big_int n)const;

        big_int getSize()const;

        big_int getP()const;
        big_int getM()const;

        std::vector<Fpelem> getElems()const;

        bool operator==(const Fp &rhs)const;
        bool operator!=(const Fp &rhs)const;


    private:
        big_int _p;
};

#endif // __FP_HPP
