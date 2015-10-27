// Implementation of a GF(p) field
#ifndef __FPXELEM_HPP
#define __FPXELEM_HPP

#include <algorithm> // find_first_of
#include "types.hpp"
#include "exceptions.hpp"
#include "fieldpElem.hpp"

// class FXelem : Felem{
class FXelem{
    public:
        //TODO Revisar si realmente copia bien el vector.
        FXelem(const std::vector<Fpelem> &v): _v(v){
            if(v.size()==0)
                throw EEmptyVector("The vector used to define the element in FXelem is empty.");
            _p=_v.back().getP();
            // Check integrity of v
            for (auto &i : _v){
                if(i.getP()!=_p)
                    throw ENotCompatible();
            // TODO check
            // Remove trailing zeros
            removeTrailingZeros(_v);
        }

        FXelem & operator=(const FXelem &rhs){
            if(&rhs != this){
                if(_p != rhs._p)
                    throw ENotCompatible();
                _v = rhs._v;
            }
            return *this;
        }

        FXelem & operator=(ll rhs){
            _v = Vector<Fpelem> = {Fpelem(rhs,_p)};
            return *this;
        }

        bool operator==(const FXelem &rhs){
            return (_v == rhs._v && _p == rhs._p);
        }

        bool operator==(ll rhs){
            return (_v.size()==1 && _v[0]==Fpelem(rhs,_p));
        }

        friend bool operator==(ll lhs, const FXelem &rhs);

        bool operator!=(const FXelem &rhs){
            return !(*this == rhs);
        }

        friend bool operator!=(ll lhs, const FXelem &rhs);

        FXelem & operator+=(const FXelem &rhs){
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
            return *this;
        }

        const FXelem operator+(const FXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return FXelem(*this) += rhs;
        }

        const FXelem operator-() const{
            std::vector<Fpelem> ret(_v.size());
            for(int i=0;i<_v.size();++i)
                ret[i]=Fpelem(-_v[i], _p);
            return FXelem(ret, _p);
        }

        FXelem & operator-=(const FXelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return (*this +=(-rhs));
        }

        const FXelem operator-(const FXelem &rhs) const{
            return FXelem(*this) -= rhs;
        }


        FXelem & operator*=(const FXelem &rhs){
            checkInSameField(rhs);

            std::vector<Fpelem> ret(rhs._v.size()+_v.size()-1,Fpelem(0,_p));
            for(int i=0;i<_v.size();++i)
                for(int j=0;j<rhs._v.size();++j){
                    ret[i+j]+=_v[i]*rhs._v[j];
                }
            _v = ret;
            return *this;
        }

        const FXelem operator*(const FXelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the *= operator
            return FXelem(*this) *= rhs;
        }

        // TODO test hard
        FXelem & operator/=(const FXelem &rhs){
            checkInSameField(rhs);
            if(rhs._v.size()==1 && rhs._v[0] == 0)
                throw EOperationUnsupported("Error. Cannot divide by the polinomial 0");
            FXelem quot(Vector(_v.size()-rhs._v.size(),0));
            unsigned int divDeg = rhs._v.size();
            Fpelem divlc = rhs._v.back();
            // Iterator to the leading coeficient
            auto remlc = _v.end();
            remlc--;
            // While the degree of the leading coefficient is greater
            //  or equal to the degree of the divisor
            while(std::distance(_v.begin(),remlc) >= divDeg){
                std::vector<Fpelem> paddingZeros(*this._v.size() - divDeg, Fpelem(0,_p));
                remlc = std::find_if(_v.rbegin(), _v.rend(),
                        std::bind1st(std::not_equal_to<Fpelem>(), Fpelem(0, _p))).base();
                paddingZeros.push_back(remlc/divlc);
                FXelem monDiv (paddingZeros);
                quot += monDiv;
                _v -= monDiv*rhs;
            }
            removeTrailingZeros(_v);

            return *this;
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
        void removeTrailingZeros(vector<Fpelem> &v){
            v.erase(
                std::find_if(v.rbegin(), v.rend(),
                    std::bind1st(std::not_equal_to<Fpelem>(), Fpelem(0, v[0]._p))).base(),
                v.end());
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

bool operator==(ll lhs, const FXelem &rhs){
    return (rhs == lhs);
}

bool operator!=(ll lhs, const FXelem &rhs){
    return rhs != lhs;
}

#endif // __FPXELEM_HPP
