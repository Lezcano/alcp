#include "fp.hpp"
#include "fpelem.hpp"
#include "fpxelem.hpp"
// types.hpp included in fp.hpp
#include "fpelem.hpp"
#include "exceptions.hpp"
#include "generalPurpose.hpp" // Miller Rabin

#include <vector>
#include <iosfwd>            // ostream
#include <string>            // to_string
#include <memory>           // unique_ptr


////////////////////////////// Fp //////////////////////////////

Fp::Fp(ll p): _p(p){
    // TODO create a look-up table for p < 2^16
    if(p<=0 || !millerRabin(p))
        throw EpNotPrime("Could not create F" + std::to_string(p) + ". " + std::to_string(p) + " is not prime.");
}

Fpelem Fp::get(ll n)const{
    return Fpelem(n, std::unique_ptr<Fp>(new Fp(*this)));
}

ll Fp::getSize()const{return _p;}

bool Fp::operator==(const Fp &rhs)const{return _p == rhs._p;}
bool Fp::operator!=(const Fp &rhs)const{return _p != rhs._p;}

std::vector<Fpelem> Fp::getElems()const{
    std::vector<Fpelem> ret;
    for(ll i=0;i<_p;++i)
        ret.push_back(this->get(i));
    return ret;
}

//////////////////////////////   Fp   //////////////////////////////

////////////////////////////// Fpelem //////////////////////////////




Fpelem::Fpelem ( const Fpelem & other){
    _num = other._num;
    _f = std::unique_ptr<Fp>(new Fp(*other._f));
}

Fpelem & Fpelem::operator=(const Fpelem &rhs){
    if(&rhs != this){
        if(*_f != *rhs._f)
            throw ENotCompatible("Fpelem assignation failed. The elements " + to_string(*this) + " and " + to_string(rhs) + " are in the fields F" + std::to_string(this->getSize()) + " and F" + std::to_string(rhs.getSize()) + " respectively.");
        _num = rhs._num;
    }
    return *this;
}

Fpelem & Fpelem::operator=(ll rhs){
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
    _num = (ll) ((_num + (ull)rhs._num)%getSize());
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

/** Russian peasant multiplication
 *   Overview
 *    It multiplies two positive integers of 63 bits and reduces them
 *     modulo p, using integers not bigger than 64 bits.
 *    It circumvents the problem of not having integers greater
 *     than 64 bits in C++.
 *    It does this by computing the multiplication adding 2^i*b(mod p)
 *     to the result if the i-th bit of a is one.
 */
Fpelem & Fpelem::operator*=(const Fpelem &rhs){
    checkInSameField(rhs);

    ull a = _num, b = rhs._num;
    ull res = 0;
    while (a != 0) {
        if (a & 1) res = (res + b) % getSize();
        a >>= 1;
        b = (b << 1) % getSize();
    }
    _num = (ll)res;
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
    ll res, aux;
    eea(_num, getSize(), res, aux);
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

int deg(const Fpelem &e){return e._num;}

// In a field the division has reminder zero
const Fpelem Fpelem::operator%(const Fpelem &rhs) const{return _f->get(0);}

ll Fpelem::getSize()const{return _f->getSize();}

const Fp Fpelem::getField()const{return *_f;}

std::string to_string(const Fpelem &e){return std::to_string(e._num);}
const Fpelem getZero(const Fpelem &e){return e.getField().get(0);}
const Fpelem getOne(const Fpelem &e){return e.getField().get(1);}

bool compatible(const Fpelem &lhs, const Fpelem &rhs){
    return lhs.getField()==rhs.getField();
}

Fpelem::Fpelem(ll num, std::unique_ptr<Fp> f): _num(num){
    _f = std::move(f);
    ll p = _f->getSize();
    _num %= p;
    if(_num < 0)
        _num += p;
}

void Fpelem::checkInSameField(const Fpelem &rhs) const{
    if(this->getField() != rhs.getField())
        throw EOperationUnsupported(
            "Error. Is not possible to add the number " + std::to_string(_num) +
            " in F" + std::to_string(getSize()) +
            " with the number " + std::to_string(rhs._num) +
            " in F" + std::to_string(rhs.getSize()));
}

Fpelem & operator+=(Fpelem &lhs, ll rhs){
    lhs+=lhs.getField().get(rhs);
    return lhs;
}

const Fpelem operator+(const Fpelem &lhs, ll rhs){
    return lhs + lhs.getField().get(rhs);
}

const Fpelem operator+(ll lhs, const Fpelem & rhs){
    return rhs.getField().get(lhs) + rhs;
}

Fpelem & operator-=(Fpelem &lhs, ll rhs){
    lhs-=lhs.getField().get(rhs);
    return lhs;
}

const Fpelem operator-(const Fpelem &lhs, ll rhs){
    return lhs - lhs.getField().get(rhs);
}

const Fpelem operator-(ll lhs, const Fpelem & rhs){
    return rhs.getField().get(lhs) - rhs;
}

bool operator==(const Fpelem & lhs, ll rhs){
    return lhs == lhs.getField().get(rhs);
}

bool operator==(ll lhs, const Fpelem &rhs){
    return rhs == lhs;
}

bool operator!=(const Fpelem & lhs, ll rhs){
    return !(lhs == rhs);
}

bool operator!=(ll lhs, const Fpelem &rhs){
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const Fpelem &e){
    os << to_string(e);
    return os;
}

////////////////////////////// Fpelem //////////////////////////////

////////////////////////////// Fxpelem //////////////////////////////


Fpxelem::Fpxelem(const Fpelem & e) : PolinomialRing<Fpxelem, Fpelem>(e){}
Fpxelem::Fpxelem(const std::vector<Fpelem> & v) : PolinomialRing<Fpxelem, Fpelem>(v){}

const Fpxelem::F Fpxelem::getField()const{
    return this->lc().getField();
}

ll Fpxelem::getSize()const{
    return this->getField().getSize();
}

// TODO comentar!!
// TODO probar!!
bool Fpxelem::irreducible()const{
    Fpxelem x({this->getField().get(0), this->getField().get(1)});
    Fpxelem xpk = x; // x^(p^k)

    for(int i=0;i<this->deg()/2;++i){
        xpk = fastPowMod<Fpxelem>(xpk, this->getSize(), *this);
        if(gcd(*this, xpk-x).deg()!=0)
            return false;
    }
    return true;
}

Fpxelem getZero(const Fpxelem &e){return Fpxelem(e.getField().get(0));}
Fpxelem getOne(const Fpxelem &e){return Fpxelem(e.getField().get(1));}


const Fpelem unit(const Fpxelem &e){ return e.lc(); }

bool compatible(const Fpxelem &lhs, const Fpxelem &rhs){
    return lhs.getField()==rhs.getField();
}

bool operator==(const Fpxelem &lhs, ll rhs){
    return lhs.deg()==0 && lhs.lc()==lhs.getField().get(rhs);
}

bool operator==(ll lhs, const Fpxelem &rhs){
    return rhs == lhs;
}

bool operator!=(const Fpxelem &lhs, ll rhs){
    return !(lhs == rhs);
}

bool operator!=(ll lhs, const Fpxelem &rhs ){
    return !(rhs == lhs);
}

////////////////////////////// Fxpelem //////////////////////////////
