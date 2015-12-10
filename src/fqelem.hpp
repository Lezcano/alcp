#ifndef __FQELEM_HPP
#define __FQELEM_HPP

#include "fq.hpp"
#include "types.hpp"

#include <iosfwd>           // ostream
#include <memory>           // unique_ptr

class Fq;

class Fqelem{
    public:
        // Base field
        using F = Fq;

        Fqelem ();
        Fqelem ( const Fqelem & );

        Fqelem & operator=(const Fqelem &rhs);
        Fqelem & operator=(big_int rhs);

        bool operator==(const Fqelem &rhs)const;

        bool operator!=(const Fqelem &rhs)const;

        Fqelem & operator+=(const Fqelem &rhs);

        Fqelem operator+(const Fqelem &rhs) const;

        Fqelem operator-() const;

        Fqelem & operator-=(const Fqelem &rhs);

        Fqelem operator-(const Fqelem &rhs) const;

        Fqelem & operator*=(const Fqelem &rhs);

        Fqelem operator*(const Fqelem &rhs) const;

        /** Multiplicative inverse */
        Fqelem inv() const;

        Fqelem & operator/=(const Fqelem &rhs);

        Fqelem operator/(const Fqelem &rhs) const;

        big_int getSize()const;

        const F getField()const;

        friend std::ostream& operator<<(std::ostream& os, const Fqelem &e);
        friend std::string to_string(const Fqelem &e);

    private:
        friend class Fq;

        Fqelem(Fpxelem n, Fq f);
        bool initialized() const;
        void checkInSameField(const Fqelem &rhs, std::string&& error) const;

        Fpxelem _num;
        Fpxelem _mod;
        std::unique_ptr<Fq> _f;
};

bool operator==(big_int lhs, const Fqelem &rhs);
bool operator==(const Fqelem &lhs, big_int rhs);
bool operator!=(big_int lhs, const Fqelem &rhs);
bool operator!=(const Fqelem &lhs, big_int rhs);
Fqelem & operator+=(Fqelem &lhs, big_int rhs);
Fqelem operator+(const Fqelem &lhs, big_int rhs);
Fqelem operator+(big_int lhs, const Fqelem & rhs);
Fqelem & operator-=(Fqelem &lhs, big_int rhs);
Fqelem operator-(const Fqelem &lhs, big_int rhs);
Fqelem operator-(big_int lhs, const Fqelem & rhs);
Fqelem & operator*=(Fqelem &lhs, big_int rhs);
Fqelem operator*(const Fqelem &lhs, big_int rhs);
Fqelem operator*(big_int lhs, const Fqelem & rhs);

bool compatible(const Fqelem &lhs, const Fqelem &rhs);
Fqelem getZero(const Fqelem &e);
Fqelem getOne(const Fqelem &e);
std::string to_string(const Fqelem &e);

#endif // __FQELEM_HPP
