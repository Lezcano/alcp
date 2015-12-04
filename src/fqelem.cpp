#include "exceptions.hpp"
#include "fqelem.hpp"
// types.hpp defined in fpelem.hpp
#include "fpelem.hpp"
#include "fpxelem.hpp"
#include "fq.hpp"
#include "generalPurpose.hpp" // ExtendedEuclideanAlgorithm (eea)

#include <iosfwd>            // ostream
#include <memory>           // unique_ptr

Fqelem::Fqelem(){ _f = nullptr; }

Fqelem::Fqelem ( const Fqelem & other) :
    _num(other._num),
    _mod(other._mod),
    _f (new Fq(*other._f))
    {}

Fqelem & Fqelem::operator=(const Fqelem &rhs){
    if(&rhs != this){
        if(!this->initialized())
            *this = Fqelem(rhs);
        else{
            checkInSameField(rhs, "Error in assignment.");
            _num = rhs._num;
        }
    }
    return *this;
}

Fqelem & Fqelem::operator=(big_int rhs){
    if(!this->initialized())
        throw std::runtime_error("Assignment to a non initialized Fqelem");
    *this = _f->get(rhs);
    return *this;
}

bool Fqelem::operator==(const Fqelem &rhs)const{
    return (_num == rhs._num && *_f == *(rhs._f));
}

bool Fqelem::operator!=(const Fqelem &rhs)const{
    return !(*this == rhs);
}

Fqelem & Fqelem::operator+=(const Fqelem &rhs){
    checkInSameField(rhs, "Addition or substraction error.");
    _num = (_num + rhs._num) % _mod;
    return *this;
}

const Fqelem Fqelem::operator+(const Fqelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the += operator
    return Fqelem(*this) += rhs;
}

const Fqelem Fqelem::operator-() const{
    return _f->get(-_num);
}

Fqelem & Fqelem::operator-=(const Fqelem &rhs){
    // We do not check if they are in the same field since
    // that will be done in the += operator
    return (*this +=(-rhs));
}

const Fqelem Fqelem::operator-(const Fqelem &rhs) const{
    return Fqelem(*this) -= rhs;
}

Fqelem & Fqelem::operator*=(const Fqelem &rhs){
    checkInSameField(rhs, "Multiplication or division error.");
    _num = (_num * rhs._num) % _mod;
    return *this;
}

const Fqelem Fqelem::operator*(const Fqelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the *= operator
    return Fqelem(*this) *= rhs;
}

/** Multiplicative inverse */
const Fqelem Fqelem::inv() const{
    if(_num == 0)
        throw EOperationUnsupported("Error. Zero has no inverse.");
    Fpxelem res = _f->get(0)._num, aux = _f->get(0)._num;
    eea(_num, _mod, res, aux);
    return _f->get(res);
}

Fqelem & Fqelem::operator/=(const Fqelem &rhs){
    // We do not check if they are in the same field since
    // that will be done in the *= operator

    return *this *= rhs.inv();
}

const Fqelem Fqelem::operator/(const Fqelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the /= operator
    return Fqelem(*this) /= rhs;
}

int deg(const Fqelem &e){return e._num.deg();}

// In a field the division has reminder zero
const Fqelem Fqelem::operator%(const Fqelem &rhs) const{return _f->get(0);}

big_int Fqelem::getSize()const{return _f->getSize();}

const Fqelem::F Fqelem::getField()const{return *_f;}

std::string to_string(const Fqelem &e){return to_string(e._num);}
const Fqelem getZero(const Fqelem &e){return e.getField().get(0);}
const Fqelem getOne(const Fqelem &e){return e.getField().get(1);}

bool compatible(const Fqelem &lhs, const Fqelem &rhs){
    return lhs.getField()==rhs.getField();
}

Fqelem::Fqelem(Fpxelem num, Fpxelem mod,  const Fq& f): _num(num), _mod(mod), _f(new Fq(f)){
    num %= mod;
}

bool Fqelem::initialized()const{ return _f != nullptr; }

void Fqelem::checkInSameField(const Fqelem &rhs, std::string&& error) const{
    if(this->getField() != rhs.getField())
        throw EOperationUnsupported(
            error + "\nThe values that caused it were " + to_string(_num) +
            " in " + to_string(this->getField()) +
            " and " + to_string(rhs._num) +
            " in F" + to_string(rhs.getField()));
}

Fqelem & operator+=(Fqelem &lhs, big_int rhs){
    lhs+=lhs.getField().get(rhs);
    return lhs;
}

const Fqelem operator+(const Fqelem &lhs, big_int rhs){
    return lhs + lhs.getField().get(rhs);
}

const Fqelem operator+(big_int lhs, const Fqelem & rhs){
    return rhs.getField().get(lhs) + rhs;
}

Fqelem & operator-=(Fqelem &lhs, big_int rhs){
    lhs-=lhs.getField().get(rhs);
    return lhs;
}

const Fqelem operator-(const Fqelem &lhs, big_int rhs){
    return lhs - lhs.getField().get(rhs);
}

const Fqelem operator-(big_int lhs, const Fqelem & rhs){
    return rhs.getField().get(lhs) - rhs;
}

Fqelem & operator*=(Fqelem &lhs, big_int rhs){
    lhs*=lhs.getField().get(rhs);
    return lhs;
}

const Fqelem operator*(const Fqelem &lhs, big_int rhs){
    return lhs * lhs.getField().get(rhs);
}

const Fqelem operator*(big_int lhs, const Fqelem & rhs){
    return rhs.getField().get(lhs) * rhs;
}

bool operator==(const Fqelem & lhs, big_int rhs){
    return lhs == lhs.getField().get(rhs);
}

bool operator==(big_int lhs, const Fqelem &rhs){
    return rhs == lhs;
}

bool operator!=(const Fqelem & lhs, big_int rhs){
    return !(lhs == rhs);
}

bool operator!=(big_int lhs, const Fqelem &rhs){
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const Fqelem &e){
    os << to_string(e);
    return os;
}
