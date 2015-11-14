// Implementation of a GF(p) field
#ifndef __POL_RING_HPP
#define __POL_RING_HPP

#include "types.hpp"
#include "exceptions.hpp"

#include <vector>
#include <algorithm>        // find_if
#include <utility>          // pair, make_pair
#include <string>           // to_string

template<typename Fxelem, typename Felem>
class PolynomialRing{
    public:
        // Inmersion from the base field
        PolynomialRing(const Felem &e): _v(std::vector<Felem>({e})){}

        PolynomialRing(const std::vector<Felem> &v): _v(v){
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
                    throw ENotCompatible("Asignation failed. The vectors "+ to_string(static_cast<Fxelem&>(*this))+ " and " + to_string(rhs) + " are not in the same ring.");
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
                        "Polynomials not in the same ring. Error when adding the polynomials " + to_string(static_cast<Fxelem&>(*this)) +
                        " and " + to_string(rhs) +  ".");
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
                        "Polynomials not in the same ring. Error when multiplying the polynomials " + to_string(static_cast<Fxelem&>(*this)) +
                        " and " + to_string(rhs) +  ".");

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
                        "Polynomials not in the same ring. Error when dividing the polynomials " + to_string(static_cast<Fxelem&>(*this)) +
                        " and " + to_string(divisor) +  ".");

            if(divisor.deg()==0 && divisor._v[0] == 0)
                throw EOperationUnsupported("Error. Cannot divide by the polynomial 0");
            if(this->deg() < divisor.deg())
                return std::make_pair(Fxelem(getZero(this->lc())), static_cast<Fxelem&>(*this));

            if(divisor.deg() == 0){
                Fxelem quot = static_cast<Fxelem&>(*this);
				Fxelem rem(static_cast<Fxelem&>(*this));
                for(auto &e : quot._v){
                    e/=divisor.lc();
                }
                return std::make_pair(quot, Fxelem(getZero(this->lc())));
            }

            Fxelem quot(getZero(this->lc()));
            Fxelem rem(static_cast<Fxelem&>(*this));

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

        // Sweet, sweet sugar
        // Adds commutativity to the operations with the base ring
        // Just with ADL one cannot add Felem + Fxelem (but you can
        //  add Fxelem + Felem...)

        friend Fxelem & operator+=(Fxelem &lhs, const Felem &rhs){
            return lhs += Fxelem(rhs);
        }

        friend const Fxelem operator+(const Fxelem &lhs, const Felem &rhs){
            return Fxelem(lhs) += rhs;
        }

        friend const Fxelem operator+(const Felem &lhs, const Fxelem &rhs){
            return Fxelem(lhs) += rhs;
        }

        friend Fxelem & operator-=(Fxelem &lhs, const Felem &rhs){
            return lhs -= Fxelem(rhs);
        }

        friend const Fxelem operator-(const Fxelem &lhs, const Felem &rhs){
            return Fxelem(lhs) -= rhs;
        }

        friend const Fxelem operator-(const Felem &lhs, const Fxelem &rhs){
            return Fxelem(lhs) -= rhs;
        }

        friend Fxelem & operator*=(Fxelem &lhs, const Felem &rhs){
            return lhs *= Fxelem(rhs);
        }

        friend const Fxelem operator*(const Fxelem &lhs, const Felem &rhs){
            return Fxelem(lhs) *= rhs;
        }

        friend const Fxelem operator*(const Felem &lhs, const Fxelem &rhs){
            return Fxelem(lhs) *= rhs;
        }

        friend Fxelem & operator/=(Fxelem &lhs, const Felem &rhs){
            return lhs /= Fxelem(rhs);
        }

        friend const Fxelem operator/(const Fxelem &lhs, const Felem &rhs){
            return Fxelem(lhs) /= rhs;
        }

        friend bool operator==(const Fxelem &lhs, Felem rhs){
            return lhs.deg()==0 && lhs.lc()==rhs;
        }

        friend bool operator==(Felem lhs, const Fxelem &rhs){
            return rhs == lhs;
        }

        friend bool operator!=(const Fxelem &lhs, Felem rhs){
            return !(lhs == rhs);
        }

        friend bool operator!=(Felem lhs, const Fxelem &rhs ){
            return !(rhs == lhs);
        }

        const Felem & operator[](int i) const {return _v[i];}
        Felem & operator[](int i) {return _v[i];}

        const Fxelem derivative()const{
           if(this->deg()==0)
               return Fxelem(getZero(this->lc()));
           std::vector<Felem> v(_v);
           for(int i=1;i<v.size();++i)
               v[i-1]=v[i]*i;
           v.pop_back();
           Fxelem ret(v);
           ret.removeTrailingZeros();

           return ret;
        }

        // Leading coefficient
        Felem lc()const{return _v.back();}

        // Degree of the polynomial
        unsigned int deg()const{return _v.size()-1;}

        // Normal form of the polynomial. It ensures the unicity of gdc for example
        friend const Fxelem normalForm(const Fxelem &e){ return e/unit(e); }

        friend std::string to_string(const Fxelem &f){
            std::string s = "";
            if(f._v.size() == 1)
                return to_string(f._v[0]);
            if(f._v.size() == 2){
                if(f._v[1] != 1)
                    s = to_string(f._v[1]);
                s += "x";
                if(f._v[0] != 0)
                    s += "+" + to_string(f._v[0]);
                return s;
            }
            if(f._v.back() != 1)
                s += to_string(f._v.back());
            s +=  "x^" + std::to_string(f._v.size()-1);

            for(int i=f._v.size()-2;i>=2;--i){
                if(f._v[i] != 0){
                    s+="+";
                    if(f._v[i] != 1)
                        s += to_string(f._v[i]);
                    s += "x^" + std::to_string(i);
                }
            }
            if(f._v[1] != 0){
                s += "+";
                if(f._v[1] != 1)
                    s += to_string(f._v[1]);
                s += "x";
            }
            if(f._v[0] != 0)
                s += "+"+ to_string(f._v[0]);
            return s;
        }

        friend std::ostream& operator<<(std::ostream& os, const Fxelem &f){
            return os << to_string(f);
        }

    protected:
        std::vector<Felem> _v;

    private:
        void removeTrailingZeros(){
            auto zero = getZero(this->lc());
            _v.erase(
                    std::find_if(
                        _v.rbegin(),
                        _v.rend(),
                        std::bind1st(std::not_equal_to<Felem>(), zero)).base(),
                    _v.end());
            // In case it was the polynomial equal to zero
            if(_v.size()==0)
                _v.push_back(zero);
        }

};



#endif //  __POL_RING_HPP
