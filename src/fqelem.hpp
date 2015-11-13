#ifndef __FQELEM_HPP
#define __FQELEM_HPP

#include "types.hpp"
// I don't think this is necessary
#include "fq.hpp"
#include <iosfwd>           // ostream
#include <memory>           // unique_ptr

class Fq;

class Fqelem{
    public:
        // Base field
        using F = Fq;

        Fqelem ( const Fqelem & );
        Fqelem & operator=(const Fqelem &rhs);

        Fqelem & operator=(big_int rhs);

        bool operator==(const Fqelem &rhs)const;

        bool operator!=(const Fqelem &rhs)const;

        Fqelem & operator+=(const Fqelem &rhs);

        const Fqelem operator+(const Fqelem &rhs) const;

        const Fqelem operator-() const;

        Fqelem & operator-=(const Fqelem &rhs);

        const Fqelem operator-(const Fqelem &rhs) const;

        Fqelem & operator*=(const Fqelem &rhs);

        const Fqelem operator*(const Fqelem &rhs) const;

        /** Multiplicative inverse */
        const Fqelem inv() const;

        Fqelem & operator/=(const Fqelem &rhs);

        const Fqelem operator/(const Fqelem &rhs) const;

        friend int deg(const Fqelem &e);

        const Fqelem operator%(const Fqelem &rhs) const;

        big_int getSize()const;

        const F getField()const;

        friend std::ostream& operator<<(std::ostream& os, const Fqelem &e);
        friend std::string to_string(const Fqelem &e);

    private:
        friend class Fq;

        Fqelem(Fpxelem n, Fpxelem mod, std::unique_ptr<F> f);
        void checkInSameField(const Fqelem &rhs) const;

        Fpxelem _num;
        std::unique_ptr<F> _f;
        Fpxelem _mod;
};

bool operator==(big_int lhs, const Fqelem &rhs);
bool operator==(const Fqelem &lhs, big_int rhs);
bool operator!=(big_int lhs, const Fqelem &rhs);
bool operator!=(const Fqelem &lhs, big_int rhs);
Fqelem & operator+=(Fqelem &lhs, big_int rhs);
const Fqelem operator+(const Fqelem &lhs, big_int rhs);
const Fqelem operator+(big_int lhs, const Fqelem & rhs);
Fqelem & operator-=(Fqelem &lhs, big_int rhs);
const Fqelem operator-(const Fqelem &lhs, big_int rhs);
const Fqelem operator-(big_int lhs, const Fqelem & rhs);

bool compatible(const Fqelem &lhs, const Fqelem &rhs);
const Fqelem getZero(const Fqelem &e);
const Fqelem getOne(const Fqelem &e);
std::string to_string(const Fqelem &e);

#endif // __FQELEM_HPP
