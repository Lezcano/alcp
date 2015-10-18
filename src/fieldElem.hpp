// Inteface of a finite field
#ifndef __FELEM_H
#define __FELEM_H

#include"types.hpp"
#include"eDElem.hpp"

class Felem : public EDelem {
    public:
        // Basic field operators
        const Felem operator+(const Felem &rhs) const = 0;
        const Felem operator-(const Felem &rhs) const = 0;
        const Felem operator*(const Felem &rhs) const = 0;
        const Felem operator/(const Felem &rhs) const = 0;
        const Felem inv() const = 0;

        // Assignment operators
        Felem & operator=(const Felem &rhs) = 0;
        Felem & operator+=(const Felem &rhs) = 0;
        Felem & operator-=(const Felem &rhs) = 0;
        Felem & operator*=(const Felem &rhs) = 0;
        Felem & operator/=(const Felem &rhs) = 0;

        // Equality operators
        bool operator==(const Felem &rhs) = 0;
        bool operator!=(const Felem &rhs) = 0;

        // Other operators
        const Felem operator-() const = 0;

        // Compatibility with ED
        const Felem operator%(const Felem &rhs) const = 0;
        int degree() const = 0;
};

#endif // __FELEM_H
