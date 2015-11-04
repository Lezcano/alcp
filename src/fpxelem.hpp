// Implementation of a GF(p) field
#ifndef __FPXELEM_HPP
#define __FPXELEM_HPP

#include <vector>
#include <algorithm>        // find_if
#include <utility>          // pair, make_pair
#include "types.hpp"
#include "exceptions.hpp"
#include "fpelem.hpp"


class FpXelem{
    public:
        // Base field
        using F = Fp;
        using Felem = Fpelem;

        FpXelem()=default; // Default ctor

        // TODO revisar move semantics
        FpXelem(const std::vector<Fpelem> &v): _v(v){
            if(v.size()==0)
                throw EEmptyVector("The vector used to define the element in FpXelem is empty.");
            _p=_v.back().getP();
            // Check integrity of v
            for (auto &i : _v)
                if(i.getP()!=_p)
                    throw ENotCompatible();
            // TODO check
            // Remove trailing zeros
            this->removeTrailingZeros();
        }

        FpXelem & operator=(const FpXelem &rhs){
            if(&rhs != this){
                if(_p != rhs._p)
                    throw ENotCompatible();
                _v = rhs._v;
            }
            return *this;
        }

        FpXelem & operator=(ll rhs){
            _v = std::vector<Fpelem>(1, Fpelem(rhs,_p));
            return *this;
        }

        bool operator==(const FpXelem &rhs)const{
            return (_v == rhs._v && _p == rhs._p);
        }

        bool operator==(ll rhs)const{
            return (_v.size()==1 && _v[0]==Fpelem(rhs,_p));
        }

        friend bool operator==(ll lhs, const FpXelem &rhs);

        bool operator!=(const FpXelem &rhs)const{
            return !(*this == rhs);
        }

        bool operator!=(ll rhs)const{
            return !(*this == rhs);
        }

        friend bool operator!=(ll lhs, const FpXelem &rhs);

        FpXelem & operator+=(const FpXelem &rhs){
            checkInSameField(rhs);
            auto v1 = _v.begin();
            auto v2 = rhs._v.begin();
            while(v1 != _v.end() && v2 != rhs._v.end()){
                *v1 += *v2;
                ++v1; ++v2;
            }
            while(v2 != rhs._v.end()){
                _v.push_back(*v2);
                ++v2;
            }
            this->removeTrailingZeros();
            return *this;
        }

        const FpXelem operator+(const FpXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return FpXelem(*this) += rhs;
        }

        const FpXelem operator-() const{
            std::vector<Fpelem> ret(_v.size());
            for(int i=0;i<_v.size();++i)
                ret[i]=Fpelem(-_v[i]);
            return FpXelem(ret);
        }

        FpXelem & operator-=(const FpXelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return (*this +=(-rhs));
        }

        const FpXelem operator-(const FpXelem &rhs) const{
            return FpXelem(*this) -= rhs;
        }


        FpXelem & operator*=(const FpXelem &rhs){
            checkInSameField(rhs);

            std::vector<Fpelem> ret(rhs._v.size()+_v.size()-1,Fpelem(0,_p));
            for(int i=0;i<_v.size();++i)
                for(int j=0;j<rhs._v.size();++j){
                    ret[i+j]+=_v[i]*rhs._v[j];
                }
            _v = ret;
            return *this;
        }

        const FpXelem operator*(const FpXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the *= operator
            return FpXelem(*this) *= rhs;
        }

        // TODO test hard
        // Implements long polynomial division
        // Return quotient and reminder in first and second respectively
        std::pair<FpXelem,FpXelem>&& div2(const FpXelem divisor){
            checkInSameField(divisor);
            if(divisor.deg()==1 && divisor._v[0] == 0)
                throw EOperationUnsupported("Error. Cannot divide by the polynomial 0");
            // Define the quotient of the corresponding size
            FpXelem quot(std::vector<Fpelem>(this->deg()-divisor.deg(),Fpelem(0,_p)));
            FpXelem rem(*this);

            // While the degree of the leading coefficient is greater
            //  or equal to the degree of the divisor
            while(rem.deg() >= divisor.deg()){
                std::vector<Fpelem> paddingZeros(rem.deg() - divisor.deg(), Fpelem(0,_p));

                paddingZeros.push_back(rem.lc()/divisor.lc());
                FpXelem monDiv (paddingZeros);
                quot += monDiv;
                rem -= monDiv*divisor;
            }
            rem.removeTrailingZeros();

            return std::move(std::make_pair(quot,rem));
        }

        FpXelem & operator/=(const FpXelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            *this = this->div2(rhs).first;
            return *this;
        }

        const FpXelem operator/(const FpXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the div2 method

            return FpXelem(*this) /= rhs;
        }

        FpXelem & operator%=(const FpXelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            *this = this->div2(rhs).second;
            return *this;
        }

        const FpXelem operator%(const FpXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the /= operator

            return FpXelem(*this) %= rhs;
        }

        friend std::ostream& operator<<(std::ostream& os, const FpXelem &f);

        // Leading coefficient
        Fpelem lc()const{return _v.back();}
        // Degree of the polynomial
        unsigned int deg()const{return _v.size()-1;}
        // Prime p of the base field Fp[X]
        ll getP()const{return _p;}

    private:
        void checkInSameField(const FpXelem &rhs) const{
            if(_p != rhs._p){
                throw EOperationUnsupported(
                    "Error when adding the polynomials " + this->to_string() +
                    " and " + rhs.to_string() +  ".");
            }
        }
        void removeTrailingZeros(){
            _v.erase(
                std::find_if(_v.rbegin(), _v.rend(),
                    std::bind1st(std::not_equal_to<Fpelem>(), Fpelem(0, _p))).base(),
                _v.end());
        }

        std::string to_string() const{
            std::string s;
            if(_v.back() != 1)
                s+= _v.back().to_string();
            s += "x^" + std::to_string(_v.size()-1);

            for(int i=_v.size()-2;i>=0;--i){
                if(_v[i]!=0){
                    s+="+";
                    if(_v[i]!=1)
                        s += _v[i].to_string();
                    s += "x^" + std::to_string(i);
                }
            }
            return s;
        }

        std::vector<Fpelem> _v;
        ll _p;
};

std::ostream& operator<<(std::ostream& os, const FpXelem &f)
{
    os << f.to_string();
    return os;
}

bool operator==(ll lhs, const FpXelem &rhs){
    return (rhs == lhs);
}

bool operator!=(ll lhs, const FpXelem &rhs){
    return (rhs != lhs);
}

typedef struct{
    FpXelem quot;
    FpXelem rem;
};

#endif // __FPXELEM_HPP
