#include <vector>
#include <algorithm>        // find_if
#include <utility>          // pair, make_pair
#include "fpxelem.hpp"
#include "types.hpp"
#include "exceptions.hpp"
#include "fpelem.hpp"
#include "fp.hpp"

// Inmersion Fp \to Fp[X]
Fpxelem::Fpxelem(const Fpelem &e): _v(std::vector<Fpelem>({e})), _f(F(e.getField())){}

// TODO delete pointer to _f!!
Fpxelem::Fpxelem(const std::vector<Fpelem> &v): _v(v), _f(F(_v.back().getField())){
    if(v.size()==0)
        throw EEmptyVector("The vector used to define the element in Fpxelem is empty.");
    ll _p=_f.getSize();
    // Check integrity of v
    for (auto &i : _v)
        if(i.getSize()!=_p)
            throw ENotCompatible("Not all the elements in the array are in the same field!");
    // Remove trailing zeros
    this->removeTrailingZeros();
}

operator std::vector< Fpxelem >()const{
    return _v;
}

Fpxelem & Fpxelem::operator=(const Fpxelem &rhs){
    if(&rhs != this){
        if(_f.getSize() != rhs._f.getSize())
            throw ENotCompatible("Asignation failed. The vectors "+ to_string()+ " and " + rhs.to_string() + " are not in the same field.");
        _v = rhs._v;
    }
    return *this;
}

bool Fpxelem::operator==(const Fpxelem &rhs)const{
    return _v == rhs._v;
}

bool Fpxelem::operator!=(const Fpxelem &rhs)const{
    return _v != rhs._v;
}


Fpxelem & Fpxelem::operator+=(const Fpxelem &rhs){
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

const Fpxelem Fpxelem::operator+(const Fpxelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the += operator
    return Fpxelem(*this) += rhs;
}

const Fpxelem Fpxelem::operator-() const{
    std::vector<Fpelem> ret(_v.size(),_f.get(0));

    for(size_t i=0;i<_v.size();++i)
        ret[i]=Fpelem(-_v[i]);
    return Fpxelem(ret);
}

Fpxelem & Fpxelem::operator-=(const Fpxelem &rhs){
    // We do not check if they are in the same field since
    // that will be done in the += operator
    return (*this +=(-rhs));
}

const Fpxelem Fpxelem::operator-(const Fpxelem &rhs) const{
    return Fpxelem(*this) -= rhs;
}

Fpxelem & Fpxelem::operator*=(const Fpxelem &rhs){
    checkInSameField(rhs);

    std::vector<Fpelem> ret(rhs._v.size()+_v.size()-1,_f.get(0));
    for(size_t i=0;i<_v.size();++i)
        for(size_t j=0;j<rhs._v.size();++j){
            ret[i+j]+=_v[i]*rhs._v[j];
        }
    _v = ret;
    return *this;
}

const Fpxelem Fpxelem::operator*(const Fpxelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the *= operator
    return Fpxelem(*this) *= rhs;
}

// TODO test hard
// Implements long polynomial division
// Return quotient and reminder in first and second respectively
std::pair<Fpxelem,Fpxelem> Fpxelem::div2(const Fpxelem &divisor){
    checkInSameField(divisor);
    if(divisor.deg()==1 && divisor._v[0] == 0)
        throw EOperationUnsupported("Error. Cannot divide by the polynomial 0");
    if(this->deg() < divisor.deg())
        return std::make_pair(Fpxelem(_f.get(0)), *this);

    // Define the quotient of the corresponding size
    Fpxelem quot(std::vector<Fpelem>(this->deg()-divisor.deg()+1,_f.get(0)));
    Fpxelem rem(*this);

    // While the degree of the leading coefficient is greater
    //  or equal to the degree of the divisor
    while(rem.deg() >= divisor.deg()){
        std::vector<Fpelem> paddingZeros(rem.deg() - divisor.deg(), _f.get(0));

        paddingZeros.push_back(rem.lc()/divisor.lc());
        Fpxelem monDiv (paddingZeros);
        quot += monDiv;
        rem -= monDiv*divisor;
    }
    rem.removeTrailingZeros();

    return std::make_pair(quot,rem);
}

Fpxelem & Fpxelem::operator/=(const Fpxelem &rhs){
    // We do not check if they are in the same field since
    // that will be done in the div2 method
    *this = this->div2(rhs).first;
    return *this;
}

const Fpxelem Fpxelem::operator/(const Fpxelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the div2 method

    return Fpxelem(*this) /= rhs;
}

Fpxelem & Fpxelem::operator%=(const Fpxelem &rhs){
    // We do not check if they are in the same field since
    // that will be done in the div2 method
    *this = this->div2(rhs).second;
    return *this;
}

const Fpxelem Fpxelem::operator%(const Fpxelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the /= operator

    return Fpxelem(*this) %= rhs;
}

const Fpelem & Fpxelem::operator[](int i) const {return _v[i];}
Fpelem & Fpxelem::operator[](int i) {return _v[i];}

// Leading coefficient
Fpelem Fpxelem::lc()const{return _v.back();}

// Degree of the polynomial
unsigned int Fpxelem::deg()const{return _v.size()-1;}

// Prime p of the base field Fp[X]
ll Fpxelem::getSize()const{return _f.getSize();}

const Fpxelem::F Fpxelem::getField()const{return _f;}

void Fpxelem::checkInSameField(const Fpxelem &rhs) const{
    if(_f != rhs._f){
        throw EOperationUnsupported(
            "Polinomials not in the same field. Error when adding the polynomials " + this->to_string() +
            " and " + rhs.to_string() +  ".");
    }
}

void Fpxelem::removeTrailingZeros(){
    _v.erase(
        std::find_if(
            _v.rbegin(),
            _v.rend(),
            std::bind1st(std::not_equal_to<Fpelem>(), _f.get(0))).base(),
        _v.end());
    // In case it was the polinomial equal to zero
    if(_v.size()==0)
        _v.push_back(_f.get(0));
}

std::string Fpxelem::to_string() const{
    std::string s;
    if(_v.size() == 1)
        return _v[0].to_string();
    if(_v.size() == 2){
        s= _v[1].to_string() + "x";
        if(_v[0] != 0)
            s+= "+" + _v[0].to_string();
        return s;
    }
    if(_v.back() != 1)
        s+=_v.back().to_string() + "x^" + std::to_string(_v.size()-1);

    for(int i=_v.size()-2;i>=2;--i){
        if(_v[i]!=0){
            s+="+";
            if(_v[i]!=1)
                s += _v[i].to_string();
            s += "x^" + std::to_string(i);
        }
    }
    if(_v[1]!=0){
        s+="+";
        if(_v[1]!=1)
            s += _v[1].to_string();
        s += "x";
    }
    if(_v[0]!=0)
        s+="+"+_v[0].to_string();
    return s;
}


std::ostream& operator<<(std::ostream& os, const Fpxelem &f){
    os << f.to_string();
    return os;
}
