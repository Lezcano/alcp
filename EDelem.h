// Interface of an arbitrary Euclidean Domain
#ifndef __EDELEM_H
#define __EDELEM_H

#include"types.h"


class EDelem{
    public:
        virtual myInt deg();
        virtual EDelem operator+(const EDelem &e1, const EDelem &e2);
        virtual EDelem operator-(const EDelem &e1, const EDelem &e2);
        virtual EDelem operator*(const EDelem &e1, const EDelem &e2);
        // Euclidean division, returns a pair (q,r) such as e1=e2*q+r
        virtual pair<EDelem,EDelem> div(const EDelem &e1, const EDelem &e2);
}

#endif // __EDELEM_H
