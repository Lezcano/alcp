#include "fpelem.hpp"
#include "fp.hpp"
#include "zelem.hpp"
#include "generalPurpose.hpp" // ExtendedEuclideanAlgorithm (eea)
#include "exceptions.hpp"

#include <stdexcept>        // std::runtime_error
#include <iosfwd>           // std::ostream
#include <string>           // std::to_string
#include <memory>           // std::unique_ptr


Fpelem::Fpelem(){ _f = nullptr; }

Fpelem::Fpelem ( const Fpelem & other) :
    _num(other._num),
    _p(other._p),
    _f(new Fp(*other._f))
    {}

Fpelem::operator big_int() const { return _num; }

Fpelem & Fpelem::operator=(const Fpelem &rhs){
    if(&rhs != this){
        if(!this->initialized())
            *this = Fpelem(rhs);
        else{
            checkInSameField(rhs, "Assignment error.");
            _num = rhs._num;
        }
    }
    return *this;
}

Fpelem & Fpelem::operator=(big_int rhs){
    if(!this->initialized())
        throw std::runtime_error("Assignment to a non initialized Fpelem");
    *this = _f->get(rhs);
    return *this;
}

bool Fpelem::operator==(const Fpelem &rhs)const{
    return (_num == rhs._num && *_f == *(rhs._f));
}

bool Fpelem::operator!=(const Fpelem &rhs)const{
    return !(*this == rhs);
}

Fpelem & Fpelem::operator+=(const Fpelem &rhs){
    checkInSameField(rhs, "Addition or substraction error.");
    this->_num = (this->_num + rhs._num) % _p;
    return *this;
}

const Fpelem Fpelem::operator+(const Fpelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the += operator
    return Fpelem(*this) += rhs;
}

const Fpelem Fpelem::operator-() const{
    return _f->get(-_num);
}

Fpelem & Fpelem::operator-=(const Fpelem &rhs){
    // We do not check if they are in the same field since
    // that will be done in the += operator
    return (*this +=(-rhs));
}

const Fpelem Fpelem::operator-(const Fpelem &rhs) const{
    return Fpelem(*this) -= rhs;
}

Fpelem & Fpelem::operator*=(const Fpelem &rhs){
    checkInSameField(rhs, "Multiplication or division error.");
    this->_num = (this->_num * rhs._num) % _p;
    return *this;
}

const Fpelem Fpelem::operator*(const Fpelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the *= operator
    return Fpelem(*this) *= rhs;
}

/** Multiplicative inverse */
const Fpelem Fpelem::inv() const{
    if(_num == 0)
        throw EOperationUnsupported("Error. Zero has no inverse.");
    big_int res, aux;
    eea(_num, _p, res, aux);
    return _f->get(res);
}

Fpelem & Fpelem::operator/=(const Fpelem &rhs){
    // We do not check if they are in the same field since
    // that will be done in the *= operator

    return *this *= rhs.inv();
}

const Fpelem Fpelem::operator/(const Fpelem &rhs) const{
    // We do not check if they are in the same field since
    // that will be done in the /= operator
    return Fpelem(*this) /= rhs;
}

big_int Fpelem::getSize()const{return this->_p;}

const Fp Fpelem::getField()const{return *_f;}

std::string to_string(const Fpelem &e){return to_string(e._num);}
const Fpelem getZero(const Fpelem &e){return e.getField().get(0);}
const Fpelem getOne(const Fpelem &e){return e.getField().get(1);}

bool compatible(const Fpelem &lhs, const Fpelem &rhs){
    return lhs.getField()==rhs.getField();
}

Fpelem::Fpelem(big_int num, const Fp& f): _num(num), _f(new Fp(f)), _p(_f->getSize()){
    _num %= _p;
    if(_num < 0)
        _num += _p;
}

bool Fpelem::initialized()const{ return _f != nullptr; }

void Fpelem::checkInSameField(const Fpelem &rhs, std::string&& error) const{
    if(this->getField() != rhs.getField())
        throw EOperationUnsupported(
            error + "\nThe values that caused it were " + to_string(_num) +
            " in F" + to_string(_p) +
            " and " + to_string(rhs._num) +
            " in F" + to_string(rhs._p));
}

Fpelem & operator+=(Fpelem &lhs, big_int rhs){
    lhs+=lhs.getField().get(rhs);
    return lhs;
}

const Fpelem operator+(const Fpelem &lhs, big_int rhs){
    return lhs + lhs.getField().get(rhs);
}

const Fpelem operator+(big_int lhs, const Fpelem & rhs){
    return rhs.getField().get(lhs) + rhs;
}

Fpelem & operator-=(Fpelem &lhs, big_int rhs){
    lhs-=lhs.getField().get(rhs);
    return lhs;
}

const Fpelem operator-(const Fpelem &lhs, big_int rhs){
    return lhs - lhs.getField().get(rhs);
}

const Fpelem operator-(big_int lhs, const Fpelem & rhs){
    return rhs.getField().get(lhs) - rhs;
}

Fpelem & operator*=(Fpelem &lhs, big_int rhs){
    lhs*=lhs.getField().get(rhs);
    return lhs;
}

const Fpelem operator*(const Fpelem &lhs, big_int rhs){
	return lhs.getField().get(rhs)*lhs ;
}

const Fpelem operator*(big_int lhs, const Fpelem & rhs){
    return rhs.getField().get(lhs) * rhs;
}

bool operator==(const Fpelem & lhs, big_int rhs){
    return lhs == lhs.getField().get(rhs);
}

bool operator==(big_int lhs, const Fpelem &rhs){
    return rhs == lhs;
}

bool operator!=(const Fpelem & lhs, big_int rhs){
    return !(lhs == rhs);
}

bool operator!=(big_int lhs, const Fpelem &rhs){
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const Fpelem &e){
    os << to_string(e);
    return os;
}
