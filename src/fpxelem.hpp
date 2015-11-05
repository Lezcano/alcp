// Implementation of a GF(p) field
#ifndef __FPXELEM_HPP
#define __FPXELEM_HPP

#include <iostream> // DEBUG
#include <vector>
#include <utility>          // pair, make_pair
#include "types.hpp"
#include "fpelem.hpp"
#include "fp.hpp"


class Fpxelem{
    public:
        // Base field
        using F = Fp;
        using Felem = Fpelem;

        Fpxelem(const Fpelem &e);
        Fpxelem(const std::vector<Fpelem> &v);

<<<<<<< HEAD
        Fpxelem & operator=(const Fpxelem &rhs);
||||||| merged common ancestors
        // TODO revisar move semantics
        // Ctor via the inmersion of Fp in Fp[X]
        Fpxelem(Fpelem e): _v({e}), _p(e.getP()){}
=======
        // TODO revisar move semantics
        // Ctor via the inmersion of Fp in Fp[X]
        Fpxelem(const Fpelem &e): _v({e}), _p(e.getP()){}

        // Default copy ctor
        Fpxelem(const Fpxelem &e)=default;
>>>>>>> af1c76ce0dadfb2bb594db1c39c8b89196283594

<<<<<<< HEAD
        bool operator==(const Fpxelem &rhs)const;
||||||| merged common ancestors
        Fpxelem(const std::vector<Fpelem> &v): _v(v){
            if(v.size()==0)
                throw EEmptyVector("The vector used to define the element in Fpxelem is empty.");
            _p=_v.back().getP();
            // Check integrity of v
            for (auto &i : _v)
                if(i.getP()!=_p)
                    throw ENotCompatible();
            // Remove trailing zeros
            this->removeTrailingZeros();
        }
=======
        Fpxelem(const std::vector<Fpelem> &v): _v(v){
            if(v.size()==0)
                throw EEmptyVector("The vector used to define the element in Fpxelem is empty.");
            _p=_v.back().getP();
            // Check integrity of v
            for (auto &i : _v)
                if(i.getP()!=_p)
                    throw ENotCompatible("Not all the elements in the array are in the same field!");
            // Remove trailing zeros
            this->removeTrailingZeros();
        }
>>>>>>> af1c76ce0dadfb2bb594db1c39c8b89196283594

<<<<<<< HEAD
        bool operator!=(const Fpxelem &rhs)const;
||||||| merged common ancestors
        Fpxelem & operator=(const Fpxelem &rhs){
            if(&rhs != this){
                if(_p != rhs._p)
                    throw ENotCompatible();
                _v = rhs._v;
            }
            return *this;
        }
=======
        Fpxelem & operator=(const Fpxelem &rhs){
            if(&rhs != this){
                if(_p != rhs._p)
                    throw ENotCompatible("Asignation failed. The vectors "+ to_string()+ " and " + rhs.to_string() + " are not in the same field.");
                _v = rhs._v;
            }
            return *this;
        }
>>>>>>> af1c76ce0dadfb2bb594db1c39c8b89196283594

        Fpxelem & operator+=(const Fpxelem &rhs);

        const Fpxelem operator+(const Fpxelem &rhs) const;

        const Fpxelem operator-() const;

        Fpxelem & operator-=(const Fpxelem &rhs);

        const Fpxelem operator-(const Fpxelem &rhs) const;

        Fpxelem & operator*=(const Fpxelem &rhs);

        const Fpxelem operator*(const Fpxelem &rhs) const;

<<<<<<< HEAD
        std::pair<Fpxelem,Fpxelem> div2(const Fpxelem &divisor);
||||||| merged common ancestors
        const Fpxelem operator-() const{
            std::vector<Fpelem> ret(_v.size());
            for(int i=0;i<_v.size();++i)
                ret[i]=Fpelem(-_v[i]);
            return Fpxelem(ret);
        }
=======
        const Fpxelem operator-() const{
            std::vector<Fpelem> ret(_v.size(),Fpelem(0,_p));

            for(size_t i=0;i<_v.size();++i)
                ret[i]=Fpelem(-_v[i]);
            return Fpxelem(ret);
        }
>>>>>>> af1c76ce0dadfb2bb594db1c39c8b89196283594

        Fpxelem & operator/=(const Fpxelem &rhs);

        const Fpxelem operator/(const Fpxelem &rhs) const;

        Fpxelem & operator%=(const Fpxelem &rhs);

        const Fpxelem operator%(const Fpxelem &rhs) const;

<<<<<<< HEAD
        const Fpelem & operator[](int i) const;
        Fpelem & operator[](int i);
||||||| merged common ancestors
            std::vector<Fpelem> ret(rhs._v.size()+_v.size()-1,Fpelem(0,_p));
            for(int i=0;i<_v.size();++i)
                for(int j=0;j<rhs._v.size();++j){
                    ret[i+j]+=_v[i]*rhs._v[j];
                }
            _v = ret;
            return *this;
        }

        const Fpxelem operator*(const Fpxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the *= operator
            return Fpxelem(*this) *= rhs;
        }

        // TODO test hard
        // Implements long polynomial division
        // Return quotient and reminder in first and second respectively
        std::pair<Fpxelem,Fpxelem>&& div2(const Fpxelem divisor){
            checkInSameField(divisor);
            if(divisor.deg()==1 && divisor._v[0] == 0)
                throw EOperationUnsupported("Error. Cannot divide by the polynomial 0");
            // Define the quotient of the corresponding size
            Fpxelem quot(std::vector<Fpelem>(this->deg()-divisor.deg(),Fpelem(0,_p)));
            Fpxelem rem(*this);

            // While the degree of the leading coefficient is greater
            //  or equal to the degree of the divisor
            while(rem.deg() >= divisor.deg()){
                std::vector<Fpelem> paddingZeros(rem.deg() - divisor.deg(), Fpelem(0,_p));

                paddingZeros.push_back(rem.lc()/divisor.lc());
                Fpxelem monDiv (paddingZeros);
                quot += monDiv;
                rem -= monDiv*divisor;
            }
            rem.removeTrailingZeros();

            return std::move(std::make_pair(quot,rem));
        }

        Fpxelem & operator/=(const Fpxelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            *this = this->div2(rhs).first;
            return *this;
        }

        const Fpxelem operator/(const Fpxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the div2 method

            return Fpxelem(*this) /= rhs;
        }

        Fpxelem & operator%=(const Fpxelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            *this = this->div2(rhs).second;
            return *this;
        }

        const Fpxelem operator%(const Fpxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the /= operator

            return Fpxelem(*this) %= rhs;
        }

        const Fpelem & operator[](int i) const {return _v[i];}
        Fpelem & operator[](int i) {return _v[i];}

        friend std::ostream& operator<<(std::ostream& os, const Fpxelem &f);
=======
            std::vector<Fpelem> ret(rhs._v.size()+_v.size()-1,Fpelem(0,_p));
            for(size_t i=0;i<_v.size();++i)
                for(size_t j=0;j<rhs._v.size();++j){
                    ret[i+j]+=_v[i]*rhs._v[j];
                }
            _v = ret;
            return *this;
        }

        const Fpxelem operator*(const Fpxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the *= operator
            return Fpxelem(*this) *= rhs;
        }

        // TODO test hard
        // Implements long polynomial division
        // Return quotient and reminder in first and second respectively
        std::pair<Fpxelem,Fpxelem> div2(const Fpxelem &divisor){
            checkInSameField(divisor);
            if(divisor.deg()==1 && divisor._v[0] == 0)
                throw EOperationUnsupported("Error. Cannot divide by the polynomial 0");
            // Define the quotient of the corresponding size
            Fpxelem quot(std::vector<Fpelem>(this->deg()-divisor.deg(),Fpelem(0,_p)));
            Fpxelem rem(*this);

            std::cout << this->deg() << " " << divisor.deg() << std::endl;
            std::cout << "quot y rem" << std::endl;
            std::cout << quot << std::endl << rem << std::endl;

            // While the degree of the leading coefficient is greater
            //  or equal to the degree of the divisor
            while(rem.deg() >= divisor.deg()){
                std::vector<Fpelem> paddingZeros(rem.deg() - divisor.deg(), Fpelem(0,_p));

                paddingZeros.push_back(rem.lc()/divisor.lc());
                Fpxelem monDiv (paddingZeros);
                quot += monDiv;
                rem -= monDiv*divisor;
            }
            rem.removeTrailingZeros();

            return std::make_pair(quot,rem);
        }

        Fpxelem & operator/=(const Fpxelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            *this = this->div2(rhs).first;
            return *this;
        }

        const Fpxelem operator/(const Fpxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the div2 method

            return Fpxelem(*this) /= rhs;
        }

        Fpxelem & operator%=(const Fpxelem &rhs){
            // We do not check if they are in the same field since
            // that will be done in the div2 method
            *this = this->div2(rhs).second;
            return *this;
        }

        const Fpxelem operator%(const Fpxelem &rhs) const{
            // We do not check if they are in the same field since
            // that will be done in the /= operator

            return Fpxelem(*this) %= rhs;
        }

        const Fpelem & operator[](int i) const {return _v[i];}
        Fpelem & operator[](int i) {return _v[i];}

        friend std::ostream& operator<<(std::ostream& os, const Fpxelem &f);
>>>>>>> af1c76ce0dadfb2bb594db1c39c8b89196283594

        // Leading coefficient
        Fpelem lc()const;
        // Degree of the polynomial
        unsigned int deg()const;
        // Prime p of the base field Fp[X]
        ll getSize()const;

        friend std::ostream& operator<<(std::ostream& os, const Fpxelem &f);

    private:
<<<<<<< HEAD
        void checkInSameField(const Fpxelem &rhs) const;
||||||| merged common ancestors
        void checkInSameField(const Fpxelem &rhs) const{
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
=======
        void checkInSameField(const Fpxelem &rhs) const{
            if(_p != rhs._p){
                throw EOperationUnsupported(
                    "Error when adding the polynomials " + this->to_string() +
                    " and " + rhs.to_string() +  ".");
            }
        }

        void removeTrailingZeros(){
            _v.erase(
                std::find_if(
                    _v.rbegin(),
                    _v.rend(),
                    std::bind1st(std::not_equal_to<Fpelem>(), Fpelem(0, _p))).base(),
                _v.end());
        }
>>>>>>> af1c76ce0dadfb2bb594db1c39c8b89196283594

        void removeTrailingZeros();

<<<<<<< HEAD
        std::string to_string() const;
||||||| merged common ancestors
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
=======
            for(int i=_v.size()-2;i>=2;--i){
                if(_v[i]!=0){
                    s+="+";
                    if(_v[i]!=1)
                        s += _v[i].to_string();
                    s += "x^" + std::to_string(i);
                }
            }
            if(_v.size()!=1){ // It cannot be zero
                if(_v[1]!=0){
                    s+="+";
                    if(_v[1]!=1)
                        s += _v[1].to_string();
                    s += "x";
                }
            }
            if(_v[0]!=0)
                s+="+"+_v[0].to_string();
            return s;
        }
>>>>>>> af1c76ce0dadfb2bb594db1c39c8b89196283594

        std::vector<Fpelem> _v;
        const F* _f;
};


#endif // __FPXELEM_HPP
