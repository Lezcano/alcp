// Implementation of a GF(p) field
#ifndef __POL_RING_HPP
#define __POL_RING_HPP

#include <vector>
#include <algorithm>        // find_if
#include <utility>          // pair, make_pair
#include <string>           // to_string

#include "types.hpp"
#include "exceptions.hpp"



template<typename Fxelem, typename Felem>
class PolinomialRing{
    public:
        // Inmersion from the base field
        PolinomialRing(const Felem &e): _v(std::vector<Felem>({e})){}

        PolinomialRing(const std::vector<Felem> &v): _v(v){
            if(v.size()==0)
                throw EEmptyVector("The vector used to define the element of R[X] is empty.");
            // Remove trailing zeros
            this->removeTrailingZeros();
            Felem aux = this->lc();
            for(auto &e : v)
                if(!compatible(aux, e))
                    throw ENotCompatible("Not all the elements in the array are in the same field.");
        }

        Fxelem & operator=(const Fxelem &rhs){
            if(&rhs != this){
                if(!compatible(this->lc(), rhs.lc()))
                    throw ENotCompatible("Asignation failed. The vectors "+ this->to_string()+ " and " + rhs.to_string() + " are not in the same ring.");
                _v = rhs._v;
            }
            return static_cast<Fxelem&>(*this);
        }

        bool operator==(const Fxelem &rhs)const{
            return _v == rhs._v;
        }

        bool operator!=(const Fxelem &rhs)const{
            return _v != rhs._v;
        }

        Fxelem & operator+=(const Fxelem &rhs){
            if(!compatible(static_cast<Fxelem&>(*this),rhs))
                throw EOperationUnsupported(
                        "Polinomials not in the same ring. Error when adding the polynomials " + this->to_string() +
                        " and " + rhs.to_string() +  ".");
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
            return static_cast<Fxelem&>(*this);
        }

        const Fxelem operator+(const Fxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return Fxelem(static_cast<const Fxelem&>(*this)) += rhs;
        }

        const Fxelem operator-() const{
            std::vector<Felem> ret(_v.size(),getZero(this->lc()));

            for(size_t i=0;i<_v.size();++i)
                ret[i]=Felem(-_v[i]);
            return Fxelem(ret);
        }

        Fxelem & operator-=(const Fxelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the += operator
            return (static_cast<Fxelem&>(*this) +=(-rhs));
        }

        const Fxelem operator-(const Fxelem &rhs) const{
            return Fxelem(static_cast<const Fxelem&>(*this)) -= rhs;
        }

        Fxelem & operator*=(const Fxelem &rhs){
            if(!compatible(static_cast<Fxelem&>(*this),rhs))
                throw EOperationUnsupported(
                        "Polinomials not in the same ring. Error when multiplying the polynomials " + this->to_string() +
                        " and " + rhs.to_string() +  ".");

            std::vector<Felem> ret(rhs._v.size()+_v.size()-1,getZero(this->lc()));
            for(size_t i=0;i<_v.size();++i)
                for(size_t j=0;j<rhs._v.size();++j){
                    ret[i+j]+=_v[i]*rhs._v[j];
                }
            _v = ret;
            return static_cast<Fxelem&>(*this);
        }

        const Fxelem operator*(const Fxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the *= operator
            return Fxelem(static_cast<const Fxelem&>(*this)) *= rhs;
        }

        // Implements long polynomial division
        // Return quotient and reminder in first and second respectively
        std::pair<Fxelem,Fxelem> div2(const Fxelem &divisor){
            if(!compatible(static_cast<Fxelem&>(*this),divisor))
                throw EOperationUnsupported(
                        "Polinomials not in the same ring. Error when dividing the polynomials " + this->to_string() +
                        " and " + divisor.to_string() +  ".");

            if(divisor.deg()==0 && divisor._v[0] == 0)
                throw EOperationUnsupported("Error. Cannot divide by the polynomial 0");
            if(this->deg() < divisor.deg())
                return std::make_pair(Fxelem(getZero(this->lc())), static_cast<Fxelem&>(*this));

            Fxelem quot(getZero(this->lc()));
            Fxelem rem(static_cast<Fxelem&>(*this));

            if(divisor.deg() == 0){
                quot = static_cast<Fxelem&>(*this);
                for(auto &e : quot._v){
                    e/=divisor.lc();
                }
                return std::make_pair(quot, Fxelem(getZero(this->lc())));
            }

            while(rem.deg() >= divisor.deg()){
                std::vector<Felem> paddingZeros(rem.deg() - divisor.deg(), getZero(this->lc()));

                paddingZeros.push_back(rem.lc()/divisor.lc());
                Fxelem monDiv (paddingZeros);
                quot += monDiv;
                rem -= monDiv*divisor;
            }
            rem.removeTrailingZeros();

            return std::make_pair(quot,rem);
        }

        Fxelem & operator/=(const Fxelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            *this = this->div2(rhs).first;
            return static_cast<Fxelem&>(*this);
        }

        const Fxelem operator/(const Fxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the div2 method

            return Fxelem(static_cast<const Fxelem&>(*this)) /= rhs;
        }

        Fxelem & operator%=(const Fxelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            *this = this->div2(rhs).second;
            return static_cast<Fxelem&>(*this);
        }

        const Fxelem operator%(const Fxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the /= operator

            return Fxelem(static_cast<const Fxelem&>(*this)) %= rhs;
        }

        const Felem & operator[](int i) const {return _v[i];}
        Felem & operator[](int i) {return _v[i];}

        // Leading coefficient
        Felem lc()const{return _v.back();}

        // Degree of the polynomial
        unsigned int deg()const{return _v.size()-1;}

        // Normal form of the polinomial. It ensures the unicity of gdc for example
        friend const Fxelem normalForm(const Fxelem &e){ return e/unit(e); }

        std::string to_string() const{
            std::string s = "";
            if(_v.size() == 1)
                return ::to_string(_v[0]);
            if(_v.size() == 2){
                if(_v[1] != 1)
                    s = ::to_string(_v[1]);
                s += "x";
                if(_v[0] != 0)
                    s += "+" + ::to_string(_v[0]);
                return s;
            }
            if(_v.back() != 1)
                s += ::to_string(_v.back());
            s +=  "x^" + std::to_string(_v.size()-1);

            for(int i=_v.size()-2;i>=2;--i){
                if(_v[i] != 0){
                    s+="+";
                    if(_v[i] != 1)
                        s += ::to_string(_v[i]);
                    s += "x^" + std::to_string(i);
                }
            }
            if(_v[1] != 0){
                s += "+";
                if(_v[1] != 1)
                    s += ::to_string(_v[1]);
                s += "x";
            }
            if(_v[0] != 0)
                s += "+"+ ::to_string(_v[0]);
            return s;
        }

        friend std::ostream& operator<<(std::ostream& os, const Fxelem &f){
            return os << f.to_string();
        }

    private:
        void removeTrailingZeros(){
            auto zero = getZero(this->lc());
            _v.erase(
                    std::find_if(
                        _v.rbegin(),
                        _v.rend(),
                        std::bind1st(std::not_equal_to<Felem>(), zero)).base(),
                    _v.end());
            // In case it was the polinomial equal to zero
            if(_v.size()==0)
                _v.push_back(zero);
        }


        std::vector<Felem> _v;
};


#endif //  __POL_RING_HPP
