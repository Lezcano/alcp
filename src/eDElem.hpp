// Interface of an arbitrary Euclidean Domain
#ifndef __EDELEM_HPP
#define __EDELEM_HPP


class EDelem {
    public:
        // Euclidean Domain Operations
        const EDelem operator+(const EDelem &rhs) const = 0;
        const EDelem operator-(const EDelem &rhs) const = 0;
        const EDelem operator*(const EDelem &rhs) const = 0;
        const EDelem operator/(const EDelem &rhs) const = 0;
        const EDelem operator%(const EDelem &rhs) const = 0;
        int degree() const = 0;
};

#endif // __EDELEM_HPP
