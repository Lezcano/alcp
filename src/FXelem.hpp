// Implementation of a GF(p) field
#ifndef __FPXELEM_HPP
#define __FPXELEM_HPP

#include "types.hpp"
#include "exceptions.hpp"
#include "fieldpElem.hpp"
//#include"fieldElem.hpp" // TODO

// class FXelem : Felem{
class FXelem{
    public:
        //TODO Revisar si realmente copia bien el vector.
        FXelem(const std::vector<Fpelem> &v): _v(v){
            /* Since the Fpelem are already checked to be prime, this is not necessary
            if(_p<=0 || !millerRabin(_p))
                throw EpNotPrime();
            */
            if(v.size()==0)
                throw EEmptyVector("The vector used to define the element in FXelem is empty.");
            // Check integrity of v
            for (auto i=v.begin(); i!=v.end();i++){
                if(i->getP()!=_p)
                    throw ENotCompatible();
            // Remove trailing zeros
            // TODO check
            value.erase(
                std::find_if(v.rbegin(), v.rend(), std::bind1st(std::not_equal_to<Fpelem>(), Fpelem(0, _v.back().p()))).base(),
                value.end());
            }
            _p=_v.back().getP();
        }

        FXelem & operator=(const FXelem &rhs){
            if(&rhs != this){
                if(_p != rhs._p)
                    throw ENotCompatible();
                _v = rhs._v;
            }
            return *this;
        }

        bool operator==(const FXelem &rhs){
            return (_v == rhs._v && _p == rhs._p);
        }

        bool operator!=(const FXelem &rhs){
            return !(*this == rhs);
        }

        FXelem & operator+=(const FXelem &rhs){
            checkInSameField(rhs);
            _num = (ll) ((_num + (ull)rhs._num)%_p);
            return *this;
        }

        const FXelem operator+(const FXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return FXelem(*this) += rhs;
        }

        const FXelem operator-() const{
            return FXelem(_p-_num, _p);
        }

        FXelem & operator-=(const FXelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return (*this +=(-rhs));
        }

        const FXelem operator-(const FXelem &rhs) const{
            return FXelem(*this) += rhs;
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
        FXelem & operator*=(const FXelem &rhs){
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

        const FXelem operator*(const FXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the *= operator
            return FXelem(*this) *= rhs;
        }

        FXelem & operator/=(const FXelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the *= operator

            return *this *= rhs.inv();
        }

        const FXelem operator/(const FXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the /= operator
            return FXelem(*this) /= rhs;
        }

        int fieldSize(){return _num;}
        const FXelem operator%(const FXelem &rhs) const{return FXelem(0,_p);}

    private:
        void checkInSameField(const FXelem &rhs) const{
            if(_p != rhs._p){
                throw EOperationUnsupported(
                    "Error when adding the polinomials " + this->to_string() +
                    " and " + rhs.to_string() "."));
            }
        }

        std::string to_string() const{
            std::string s;
            if(_v.back() == 1)
                s += "x^" + std::to_string(_v.size());
            else if(_v.back() == -1)
                s += "-x^" + std::to_string(_v.size());
            else
                s += std::to_string(_v.back()) + "-x^" + std::to_string(_v.size());

            for(int i=_v.size()-2;i>=0;--i){
                if(_v[i]!=0){
                    if(_v[i] > 0)
                        s+= "+";
                    s+=std::to_string(_v[i]) + "x^" + std::to_string(i);
                }
            }
            s += "(mod " + std::to_string(_p) + ")";
            return s;
        }

        vector<fieldpElem> _v;
        ll _p;
};

#endif // __FPXELEM_HPP
