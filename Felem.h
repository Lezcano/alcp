// Inteface of a finite field
#ifndef __FELEM_H
#define __FELEM_H

#include"types.h"

class Felem : public EDelem{
    public:
        virtual Felem operator+(const Felem &f1, const Felem &f2) const = 0;
        virtual Felem operator-(const Felem &f1, const Felem &f2) const = 0;
        virtual Felem operator*(const Felem &f1, const Felem &f2) const = 0;
        virtual Felem operator/(const Felem &f1, const Felem &f2) const = 0;
}

#endif // __FELEM_H
