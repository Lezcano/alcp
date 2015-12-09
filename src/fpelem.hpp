#ifndef __FPELEM_HPP
#define __FPELEM_HPP

#include "fp.hpp"
#include "types.hpp"

#include <iosfwd>           // ostream
#include <memory>           // unique_ptr

class Fp;

class Fpelem{
    public:
        // Base field
        using F = Fp;

        Fpelem ();
        Fpelem ( const Fpelem & );
        explicit operator big_int() const;//TODO

        Fpelem & operator=(const Fpelem &rhs);
        Fpelem & operator=(big_int rhs);

        bool operator==(const Fpelem &rhs)const;

        bool operator!=(const Fpelem &rhs)const;

        Fpelem & operator+=(const Fpelem &rhs);

        Fpelem operator+(const Fpelem &rhs) const;

        Fpelem operator-() const;

        Fpelem & operator-=(const Fpelem &rhs);

        Fpelem operator-(const Fpelem &rhs) const;

        Fpelem & operator*=(const Fpelem &rhs);

        Fpelem operator*(const Fpelem &rhs) const;

        /** Multiplicative inverse */
        Fpelem inv() const;

        Fpelem & operator/=(const Fpelem &rhs);

        Fpelem operator/(const Fpelem &rhs) const;

        big_int getSize()const;

        const F getField()const;

        friend std::ostream& operator<<(std::ostream& os, const Fpelem &e);
        friend std::string to_string(const Fpelem &e);

    private:
        friend class Fp;

        Fpelem(big_int num, Fp f);
        bool initialized() const;
        void checkInSameField(const Fpelem &rhs, std::string&& error) const;

        big_int _num;
        big_int _p;
        std::unique_ptr<F> _f;
};

bool operator==(big_int lhs, const Fpelem &rhs);
bool operator==(const Fpelem &lhs, big_int rhs);
bool operator!=(big_int lhs, const Fpelem &rhs);
bool operator!=(const Fpelem &lhs, big_int rhs);
Fpelem & operator+=(Fpelem &lhs, big_int rhs);
Fpelem operator+(const Fpelem &lhs, big_int rhs);
Fpelem operator+(big_int lhs, const Fpelem & rhs);
Fpelem & operator-=(Fpelem &lhs, big_int rhs);
Fpelem operator-(const Fpelem &lhs, big_int rhs);
Fpelem operator-(big_int lhs, const Fpelem & rhs);
Fpelem & operator*=(Fpelem &lhs, big_int rhs);
Fpelem operator*(const Fpelem &lhs, big_int rhs);
Fpelem operator*(big_int lhs, const Fpelem & rhs);

bool compatible(const Fpelem &lhs, const Fpelem &rhs);
Fpelem getZero(const Fpelem &e);
Fpelem getOne(const Fpelem &e);
std::string to_string(const Fpelem &e);

#endif // __FPELEM_HPP
