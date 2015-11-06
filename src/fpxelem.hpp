// Implementation of a GF(p) field
#ifndef __FPXELEM_HPP
#define __FPXELEM_HPP

#include <vector>
#include <utility>          // pair, make_pair
#include "types.hpp"
#include "fpelem.hpp"
#include "fp.hpp"


class Fpxelem{
    public:
        // Base field
        using F = Fp;
        using Felem = Fpelem;

        Fpxelem(const Fpelem &e);
        Fpxelem(const std::vector<Fpelem> &v);

        operator std::vector< Fpxelem >()const;

        Fpxelem & operator=(const Fpxelem &rhs);

        bool operator==(const Fpxelem &rhs)const;
        bool operator==(ll rhs)const;

        bool operator!=(const Fpxelem &rhs)const;
        bool operator!=(ll rhs)const;

        Fpxelem & operator+=(const Fpxelem &rhs);

        const Fpxelem operator+(const Fpxelem &rhs) const;

        const Fpxelem operator-() const;

        Fpxelem & operator-=(const Fpxelem &rhs);

        const Fpxelem operator-(const Fpxelem &rhs) const;

        Fpxelem & operator*=(const Fpxelem &rhs);

        const Fpxelem operator*(const Fpxelem &rhs) const;

        std::pair<Fpxelem,Fpxelem> div2(const Fpxelem &divisor);

        Fpxelem & operator/=(const Fpxelem &rhs);

        const Fpxelem operator/(const Fpxelem &rhs) const;

        Fpxelem & operator%=(const Fpxelem &rhs);

        const Fpxelem operator%(const Fpxelem &rhs) const;

        const Fpelem & operator[](int i) const;
        Fpelem & operator[](int i);

        // Leading coefficient
        Fpelem lc()const;
        // Degree of the polynomial
        unsigned int deg()const;
        // Prime p of the base field Fp[X]
        ll getSize()const;

        const F getField()const;

        friend std::ostream& operator<<(std::ostream& os, const Fpxelem &f);

    private:
        void checkInSameField(const Fpxelem &rhs) const;

        void removeTrailingZeros();

        std::string to_string() const;

        std::vector<Fpelem> _v;
        const F _f;
};

const Fpxelem unit(const Fpxelem &f);
const Fpxelem normalForm(const Fpxelem &f);


#endif // __FPXELEM_HPP
