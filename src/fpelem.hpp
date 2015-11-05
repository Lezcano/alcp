#ifndef __FPELEM_HPP
#define __FPELEM_HPP

#include<iosfwd>            // ostream
#include "types.hpp"
#include "fp.hpp"

class Fp;

class Fpelem{
    public:
<<<<<<< HEAD
        // Base field
        typedef Fp F;
||||||| merged common ancestors
        Fpelem()=default; // Default ctor
        Fpelem(ll num, ll p){
            if(p<=0 || !millerRabin(p))
                throw EpNotPrime();
            _p = p;
            _num = num % p;
            if(_num < 0)
                _num = (_num + p) % p;
        }

        Fpelem & operator=(const Fpelem &rhs){
            if(&rhs != this){
                if(_p != rhs._p)
                    throw ENotCompatible();
                _num = rhs._num;
            }
            return *this;
        }

        bool operator==(const Fpelem &rhs)const{
            return (_num == rhs._num && _p == rhs._p);
        }

        friend bool operator==(ll lhs, const Fpelem &rhs);
        friend bool operator==(const Fpelem &lhs, ll rhs);

        bool operator!=(const Fpelem &rhs)const{
            return !(*this == rhs);
        }

        friend bool operator!=(ll lhs, const Fpelem &rhs);
        friend bool operator!=(const Fpelem &lhs, ll rhs);

        Fpelem & operator+=(const Fpelem &rhs){
            checkInSameField(rhs);
            _num = (ll) ((_num + (ull)rhs._num)%_p);
            return *this;
        }

        const Fpelem operator+(const Fpelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return Fpelem(*this) += rhs;
        }
        friend const Fpelem operator+(const Fpelem &lhs, ll rhs);
        friend const Fpelem operator+(ll lhs, const Fpelem & rhs);

        const Fpelem operator-() const{
            return Fpelem(_p-_num, _p);
        }

        Fpelem & operator-=(const Fpelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return (*this +=(-rhs));
        }

        const Fpelem operator-(const Fpelem &rhs) const{
            return Fpelem(*this) -= rhs;
        }
        friend const Fpelem operator-(const Fpelem &lhs, ll rhs);
        friend const Fpelem operator-(ll lhs, const Fpelem & rhs);


        /** Russian peasant multiplication
         *   Overview
         *    It multiplies two positive integers of 63 bits and reduces them
         *     modulo p, using integers not bigger than 64 bits.
         *    It circunvents the problem of not having integers greater
         *     than 64 bits in C++.
         *    It does this by computing the multiplication adding 2^i*b(mod p)
         *     to the result if the i-th bit of a is one.
         */
        Fpelem & operator*=(const Fpelem &rhs){
            checkInSameField(rhs);

            ull a = _num, b = rhs._num;
            ull res = 0;
            while (a != 0) {
                if (a & 1) res = (res + b) % _p;
                a >>= 1;
                b = (b << 1) % _p;
            }
            _num = (ll)res;
            return *this;
        }

        const Fpelem operator*(const Fpelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the *= operator
            return Fpelem(*this) *= rhs;
        }
=======
        Fpelem()=default;               // Default ctor
        Fpelem(const Fpelem&)=default;  // Default copy ctor
        Fpelem(ll num, ll p){
            // TODO create a look-up list for p < 2^16
            if(p<=0 || !millerRabin(p))
                throw EpNotPrime("Could not create Fpelem. The parameter p is not prime.");
            _p = p;
            _num = num % p;
            if(_num < 0)
                _num = (_num + p) % p;
        }

        Fpelem & operator=(const Fpelem &rhs){
            if(&rhs != this){
                if(_p != rhs._p)
                    throw ENotCompatible("Fpelem assignation failed. The elements " + this->to_string() + " and " + rhs.to_string()
                                         + " are in the fields F" + std::to_string(_p) + " and F" + std::to_string(rhs._p) + " respectively.");
                _num = rhs._num;
            }
            return *this;
        }

        bool operator==(const Fpelem &rhs)const{
            return (_num == rhs._num && _p == rhs._p);
        }

        friend bool operator==(ll lhs, const Fpelem &rhs);
        friend bool operator==(const Fpelem &lhs, ll rhs);

        bool operator!=(const Fpelem &rhs)const{
            return !(*this == rhs);
        }

        friend bool operator!=(ll lhs, const Fpelem &rhs);
        friend bool operator!=(const Fpelem &lhs, ll rhs);

        Fpelem & operator+=(const Fpelem &rhs){
            checkInSameField(rhs);
            _num = (ll) ((_num + (ull)rhs._num)%_p);
            return *this;
        }

        const Fpelem operator+(const Fpelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return Fpelem(*this) += rhs;
        }
        friend const Fpelem operator+(const Fpelem &lhs, ll rhs);
        friend const Fpelem operator+(ll lhs, const Fpelem & rhs);

        const Fpelem operator-() const{
            return Fpelem(_p-_num, _p);
        }

        Fpelem & operator-=(const Fpelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return (*this +=(-rhs));
        }

        const Fpelem operator-(const Fpelem &rhs) const{
            return Fpelem(*this) -= rhs;
        }
        friend const Fpelem operator-(const Fpelem &lhs, ll rhs);
        friend const Fpelem operator-(ll lhs, const Fpelem & rhs);


        /** Russian peasant multiplication
         *   Overview
         *    It multiplies two positive integers of 63 bits and reduces them
         *     modulo p, using integers not bigger than 64 bits.
         *    It circunvents the problem of not having integers greater
         *     than 64 bits in C++.
         *    It does this by computing the multiplication adding 2^i*b(mod p)
         *     to the result if the i-th bit of a is one.
         */
        Fpelem & operator*=(const Fpelem &rhs){
            checkInSameField(rhs);

            ull a = _num, b = rhs._num;
            ull res = 0;
            while (a != 0) {
                if (a & 1) res = (res + b) % _p;
                a >>= 1;
                b = (b << 1) % _p;
            }
            _num = (ll)res;
            return *this;
        }

        const Fpelem operator*(const Fpelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the *= operator
            return Fpelem(*this) *= rhs;
        }
>>>>>>> af1c76ce0dadfb2bb594db1c39c8b89196283594

        Fpelem(const Fpelem&);


        Fpelem & operator=(const Fpelem &rhs);

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

        const F* getField()const;

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
const Fpelem operator+(const Fpelem &lhs, ll rhs);
const Fpelem operator+(ll lhs, const Fpelem & rhs);
const Fpelem operator-(const Fpelem &lhs, ll rhs);
const Fpelem operator-(ll lhs, const Fpelem & rhs);

#endif // __FPELEM_HPP
