#ifndef __FPELEM_HPP
#define __FPELEM_HPP

#include<iosfwd>            // ostream
#include "types.hpp"
#include "fp.hpp"

class Fp;

class Fpelem{
    public:
        // Base field
        typedef Fp F;
        Fpelem(const Fpelem&);

        Fpelem & operator=(const Fpelem &rhs);

        Fpelem & operator=(ll rhs);

        bool operator==(const Fpelem &rhs)const;

        bool operator!=(const Fpelem &rhs)const;

        Fpelem & operator+=(const Fpelem &rhs);

        const Fpelem operator+(const Fpelem &rhs) const;

        const Fpelem operator-() const;

        Fpelem & operator-=(const Fpelem &rhs);

        const Fpelem operator-(const Fpelem &rhs) const;

        Fpelem & operator*=(const Fpelem &rhs);

        const Fpelem operator*(const Fpelem &rhs) const;

        /** Multiplicative inverse */
        const Fpelem inv() const;

        Fpelem & operator/=(const Fpelem &rhs);

        const Fpelem operator/(const Fpelem &rhs) const;

        friend int deg(const Fpelem &e);

        const Fpelem operator%(const Fpelem &rhs) const;

        ll getSize()const;

        const F getField()const;

        std::string to_string()const;

        friend std::ostream& operator<<(std::ostream& os, const Fpelem &e);

    private:
        friend class Fp;

        Fpelem(ll num, const Fp* f);
        void checkInSameField(const Fpelem &rhs) const;

        ll _num;
        const F* _f;
};
bool operator==(ll lhs, const Fpelem &rhs);
bool operator==(const Fpelem &lhs, ll rhs);
bool operator!=(ll lhs, const Fpelem &rhs);
bool operator!=(const Fpelem &lhs, ll rhs);
Fpelem & operator+=(Fpelem &lhs, ll rhs);
const Fpelem operator+(const Fpelem &lhs, ll rhs);
const Fpelem operator+(ll lhs, const Fpelem & rhs);
Fpelem & operator-=(Fpelem &lhs, ll rhs);
const Fpelem operator-(const Fpelem &lhs, ll rhs);
const Fpelem operator-(ll lhs, const Fpelem & rhs);

#endif // __FPELEM_HPP
