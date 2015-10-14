// Implementation of a GF(p) field
#ifndef __FPELEM_H
#define __FPELEM_H

#include"types.h"
#include"exceptions.h"
#include"generalPurpose.h" // Miller Rabin, eea

class Fpelem : public Felem{
    public:
        Fpelem(ll num, ll p){
            if(p<=0 || !millerRabin(p))
                throw EpNotPrime();
            if(num < 0)
        }

        Fpelem& operator=(const Fpelem &rhs){
            if(&rhs != this){
                _num = rhs._num;
                _p = rhs._p;
            }
            return *this;
        }

        virtual Fpelem operator+(const Fpelem &f1, const Fpelem &f2) const;
        virtual Fpelem operator-(const Fpelem &f1, const Fpelem &f2) const;
        Fpelem operator*(const Fpelem &f1, const Fpelem &f2) const{
            if(f1._p != f2._p)
                throw EOperationUnsupported(
                    "Error when adding the number " + to_string(f1._num) +
                    " in F" + to_string(f1._p)+
                    " with the number " + to_string(f2._num) +
                    " in F" + to_string(f2._p));

        }
        virtual Fpelem operator/(const Fpelem &f1, const Fpelem &f2) const;
    private:
        ll _num;
        ll _p;
}

#endif // __FPELEM_H
