// Implementation of a GF(p) field
#ifndef __FPXELEM_HPP
#define __FPXELEM_HPP

#include <algorithm> // find_first_of
#include "types.hpp"
#include "exceptions.hpp"
#include "fieldpElem.hpp"

class FpXelem{
    public:
        //TODO Revisar si realmente copia bien el vector.
        FpXelem(const std::vector<Fpelem> &v): _v(v){
            if(v.size()==0)
                throw EEmptyVector("The vector used to define the element in FpXelem is empty.");
            _p=_v.back().getP();
            // Check integrity of v
            for (auto &i : _v){
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
            _v = Vector<Fpelem> = {Fpelem(rhs,_p)};
            return *this;
        }

        bool operator==(const FpXelem &rhs){
            return (_v == rhs._v && _p == rhs._p);
        }

        bool operator==(ll rhs){
            return (_v.size()==1 && _v[0]==Fpelem(rhs,_p));
        }

        friend bool operator==(ll lhs, const FpXelem &rhs);

        bool operator!=(const FpXelem &rhs){
            return !(*this == rhs);
        }

        friend bool operator!=(ll lhs, const FpXelem &rhs);

        FpXelem & operator+=(const FpXelem &rhs){
            checkInSameField(rhs);
            auto v1 = _v.begin(), v2 = rhs.begin();
            auto e1 = _v.end(), e2 = rhs.end();
            while(v1 != e1 && v2 != e2){
                *v1 += *v2;
                ++v1; ++v2;
            }
            while(v2 != e2){
                v1.push_back(*v2);
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
                ret[i]=Fpelem(-_v[i], _p);
            return FpXelem(ret, _p);
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
        // Implements long polinomial division
        FpXdiv&& div2(const FpXelem divisor){
            checkInSameField(divisor);
            if(divisor._v.size()==1 && divisor._v[0] == 0)
                throw EOperationUnsupported("Error. Cannot divide by the polinomial 0");
            FpXdiv ret{
                FpXelem(Vector(_v.size()-divisor._v.size(),0)), //quot
                FpXelem(*this) //rem
            };

            // While the degree of the leading coefficient is greater
            //  or equal to the degree of the divisor
            while(ret.rem.deg() >= divisor.deg()){
                std::vector<Fpelem> paddingZeros(remDeg - divDeg, Fpelem(0,_p));

                paddingZeros.push_back(ret.rem.lc()/divisor.lc());
                FpXelem monDiv (paddingZeros);
                ret.quot += monDiv;
                ret.rem -= monDiv*divisor;
            }
            ret.rem.removeTrailingZeros();

            return ret;
        }

        FpXelem & operator/=(const FpXelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            _v = div2(_v, rhs).quot;
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
            _v = div2(_v, rhs).rem;
            return *this;
        }

        const FpXelem operator%(const FpXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the /= operator

            return FpXelem(*this) %= rhs;
        }

        ostream& operator<<(ostream& os, const FpXelem &f)
        {
            os << f.to_string();
            return os;
        }

        // Leading coefficient
        fieldpElem lc(){return _v.back();}
        // Degree of the polinomial
        unsigned int deg(){return _v.size()-1;}
        // Prime p of the base field Fp[X]
        ll getP(){return _p;}

    private:
        void checkInSameField(const FpXelem &rhs) const{
            if(_p != rhs._p){
                throw EOperationUnsupported(
                    "Error when adding the polinomials " + this->to_string() +
                    " and " + rhs.to_string() "."));
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
            return s;
        }

        vector<fieldpElem> _v;
        ll _p;
};

bool FpXelem::operator==(ll lhs, const FpXelem &rhs){
    return (rhs == lhs);
}

bool FpXelem::operator!=(ll lhs, const FpXelem &rhs){
    return rhs != lhs;
}

#endif // __FPXELEM_HPP
