#ifndef __FIELDQ_HPP
#define __FIELDQ_HPP

#include "types.hpp"
#include "fp.hpp"
#include "fpxelem.hpp"

#include <iosfwd>           // ostream
#include <memory>           // unique_ptr
#include <vector>


////////////////////////////// Fq //////////////////////////////
class Fqelem;

class Fq{
    public:
        Fq(ll p, int n);

        Fqelem get(ll n)const;

        Fqelem get(Fpxelem f)const;

        ll getSize()const;

        std::vector<Fqelem> getElems()const;

        bool operator==(const Fq &rhs)const;
        bool operator!=(const Fq &rhs)const;

        friend std::string to_string(const Fq &e);
    private:

        ll _p;
        int _n;
        Fp _base;
        Fpxelem _mod;
};



////////////////////////////// Fqelem //////////////////////////////

class Fqelem{
    public:
        // Base field
        using F = Fq;

        Fqelem ( const Fqelem & );
        Fqelem & operator=(const Fqelem &rhs);

        Fqelem & operator=(ll rhs);

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

        ll getSize()const;

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

bool operator==(ll lhs, const Fqelem &rhs);
bool operator==(const Fqelem &lhs, ll rhs);
bool operator!=(ll lhs, const Fqelem &rhs);
bool operator!=(const Fqelem &lhs, ll rhs);
Fqelem & operator+=(Fqelem &lhs, ll rhs);
const Fqelem operator+(const Fqelem &lhs, ll rhs);
const Fqelem operator+(ll lhs, const Fqelem & rhs);
Fqelem & operator-=(Fqelem &lhs, ll rhs);
const Fqelem operator-(const Fqelem &lhs, ll rhs);
const Fqelem operator-(ll lhs, const Fqelem & rhs);

bool compatible(const Fqelem &lhs, const Fqelem &rhs);
const Fqelem getZero(const Fqelem &e);
const Fqelem getOne(const Fqelem &e);
std::string to_string(const Fqelem &e);

////////////////////////////// Fqelem //////////////////////////////



#endif // __FQ_HPP
