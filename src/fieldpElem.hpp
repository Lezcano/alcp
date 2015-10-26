// Implementation of a GF(p) field
#ifndef __FPELEM_HPP
#define __FPELEM_HPP

#include<iosfwd>
#include"types.hpp"
#include"exceptions.hpp"
#include"generalPurpose.hpp" // Miller Rabin, eea

// class Fpelem : Felem{
class Fpelem{
    public:
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

        bool operator==(const Fpelem &rhs){
            return (_num == rhs._num && _p == rhs._p);
        }

        bool operator!=(const Fpelem &rhs){
            return !(*this == rhs);
        }

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

        const Fpelem operator-() const{
            return Fpelem(_p-_num, _p);
        }

        Fpelem & operator-=(const Fpelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return (*this +=(-rhs));
        }

        const Fpelem operator-(const Fpelem &rhs) const{
            return Fpelem(*this) += rhs;
        }


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

        /** Multiplicative inverse */
        const Fpelem inv() const{
            if(_num == 0)
                throw EOperationUnsupported("Error. Zero has no inverse.");
            ll res, aux;
            eea(_num, _p, res, aux);
            return Fpelem(res, _p);
        }

        Fpelem & operator/=(const Fpelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the *= operator

            return *this *= rhs.inv();
        }

        const Fpelem operator/(const Fpelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the /= operator
            return Fpelem(*this) /= rhs;
        }

        int degree(){return _num;}
        const Fpelem operator%(const Fpelem &rhs) const{return Fpelem(0,_p);}

        ll getP(){return _p;}

        ostream& operator<<(ostream& os, const Fpelem &f)
        {
            os << f._num;
            return os;
        }

    private:
        void checkInSameField(const Fpelem &rhs) const{
            if(_p != rhs._p)
                throw EOperationUnsupported(
                    "Error. Is not possible to add the number " + std::to_string(_num) +
                    " in F" + std::to_string(_p) +
                    " with the number " + std::to_string(rhs._num) +
                    " in F" + std::to_string(rhs._p));
        }

        ll _num;
        ll _p;
};

#endif // __FPELEM_HPP
