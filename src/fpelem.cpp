#include "fpelem.hpp"
#include "fp.hpp"
#include "zelem.hpp"
#include "generalPurpose.hpp" // ExtendedEuclideanAlgorithm (eea)
#include "exceptions.hpp"

#include <iosfwd>            // ostream
#include <string>            // to_string
#include <memory>           // unique_ptr

Fpelem::Fpelem ( const Fpelem & other){
    _num = other._num;
    _p = other._p;
    _f = std::unique_ptr<Fp>(new Fp(*other._f));
}

Fpelem & Fpelem::operator=(const Fpelem &rhs){
    if(&rhs != this){
        if(*_f != *rhs._f)
            throw ENotCompatible("Fpelem assignation failed. The elements " + to_string(*this) + " and " + to_string(rhs) + " are in the fields F" + to_string(this->_p) + " and F" + to_string(rhs._p) + " respectively.");
        _num = rhs._num;
    }
    return *this;
}

Fpelem & Fpelem::operator=(big_int rhs){
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
    checkInSameField(rhs);
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
    checkInSameField(rhs);
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

big_int Fpelem::getSize()const{return _f->_p;}

const Fp Fpelem::getField()const{return *_f;}

std::string to_string(const Fpelem &e){return to_string(e._num);}
const Fpelem getZero(const Fpelem &e){return e.getField().get(0);}
const Fpelem getOne(const Fpelem &e){return e.getField().get(1);}

bool compatible(const Fpelem &lhs, const Fpelem &rhs){
    return lhs.getField()==rhs.getField();
}

Fpelem::Fpelem(big_int num, std::unique_ptr<Fp> f): _num(num){
    _f = std::move(f);
    _p = _f->getSize();
    _num %= _p;
    if(_num < 0)
        _num += _p;
}

void Fpelem::checkInSameField(const Fpelem &rhs) const{
    if(this->getField() != rhs.getField())
        throw EOperationUnsupported(
            "Error. Is not possible to add the number " + to_string(_num) +
            " in F" + to_string(_p) +
            " with the number " + to_string(rhs._num) +
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
	std::cout << lhs << std::endl;
	std::cout << rhs << std::endl;
	std::cout << lhs.getField().get(rhs) << std::endl;
	std::cout << lhs*lhs.getField().get(rhs)<< std::endl;
	lhs*lhs.getField().get(rhs);
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
