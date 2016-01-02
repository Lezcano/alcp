#ifndef __FELEMBASE_HPP
#define __FELEMBASE_HPP

#include <stdexcept>        // std::runtime_error
#include <iosfwd>           // std::ostream
#include <memory>           // std::unique_ptr,  std::make_unique
#include <string>           // std::to_string
#include <type_traits>      // std::enable_if, std::is_integral

#include "generalPurpose.hpp" // ExtendedEuclideanAlgorithm (eea)
#include "exceptions.hpp"

namespace alcp {
    template<template<class> class FelemBase,
            template<class> class Fbase, class Quotient, class Integer>
    class QuotientRing {
    public:
        // Base field
        using F = Fbase<Integer>;
        using Felem = FelemBase<Integer>;

        QuotientRing() : _f(nullptr) { }

        QuotientRing(const QuotientRing &other) :
                _num(other._num),
                _mod(other._mod),
                _f(new F(*other._f)) { }

        // Generalized copy constructor
        template<class Int2, class = typename std::enable_if<std::is_integral<Int2>::value>::type>
        QuotientRing(const QuotientRing<FelemBase, Fbase, Quotient, Int2> &other) :
                _num(other._num),
                _mod(other._mod),
                _f(new F(*other._f)) { }

        QuotientRing(QuotientRing &&) = default;

        explicit operator Quotient() const { return _num; }

        // Copy assignment
        QuotientRing &operator=(const QuotientRing &rhs) {
            if (&rhs != this) {
                if (!this->initialized()) {
                    _mod = rhs._mod;
                    _f.reset(new F(*rhs._f));
                }
                else {
                    checkInSameField(rhs, "Error in assignment.");
                }
                _num = rhs._num;
            }
            return static_cast<Felem &>(*this);
        }

        // Generalized copy assignment
        template<class Int2, class = typename std::enable_if<std::is_integral<Int2>::value>::type>
        QuotientRing &operator=(const QuotientRing<FelemBase, Fbase, Quotient, Int2> &rhs) {
            if (&rhs != this) {
                if (!this->initialized()) {
                    _mod = rhs._mod;
                    _f.reset(new F(*rhs._f));
                }
                else {
                    checkInSameField(rhs, "Error in assignment.");
                }
                _num = rhs._num;
            }
            return static_cast<Felem &>(*this);
        }

        // Generalized move assignment
        template<class Int2, class = typename std::enable_if<std::is_integral<Int2>::value>::type>
        QuotientRing &operator=(QuotientRing<FelemBase, Fbase, Quotient, Int2> &&rhs) {
            if (&rhs != this) {
                if (!this->initialized()) {
                    _mod = std::move(rhs._mod);
                    _f.reset(new F(*rhs._f));
                    rhs._f.reset(nullptr);
                }
                else {
                    checkInSameField(rhs, "Error in assignment.");
                }
                _num = std::move(rhs._num);
            }
            return static_cast<Felem &>(*this);
        }

        template<class Int2, class = typename std::enable_if<std::is_integral<Int2>::value>::type>
        Felem &operator=(Int2 rhs) {
            if (!this->initialized())
                throw std::runtime_error("Assignment to a non initialized Felem");
            static_cast<Felem &>(*this) = _f->get(rhs);
            return static_cast<Felem &>(*this);
        }

        bool operator==(const Felem &rhs) const {
            return (_num == rhs._num && *_f == *(rhs._f));
        }

        bool operator!=(const Felem &rhs) const {
            return !(static_cast<const Felem &>(*this) == rhs);
        }

        Felem &operator+=(const Felem &rhs) {
            checkInSameField(rhs, "Addition or substraction error.");
            _num = (_num + rhs._num) % _mod;
            return static_cast<Felem &>(*this);
        }

        Felem operator+(const Felem &rhs) const {
            return Felem(static_cast<const Felem &>(*this)) += rhs;
        }

        Felem operator-() const {
            return _f->get(-_num);
        }

        Felem &operator-=(const Felem &rhs) {
            return (static_cast<Felem &>(*this) += (-rhs));
        }

        Felem operator-(const Felem &rhs) const {
            return Felem(static_cast<const Felem &>(*this)) -= rhs;
        }

        Felem &operator*=(const Felem &rhs) {
            checkInSameField(rhs, "Multiplication or division error.");
            _num = (_num * rhs._num) % _mod;
            return static_cast<Felem &>(*this);
        }

        Felem operator*(const Felem &rhs) const {
            return Felem(static_cast<const Felem &>(*this)) *= rhs;
        }

        /** Multiplicative inverse */
        Felem inv() const {
            if (_num == 0)
                throw EOperationUnsupported("Error. Zero has no inverse.");
            Quotient res, aux;
            eea(_num, _mod, res, aux);
            return _f->get(res);
        }

        Felem &operator/=(const Felem &rhs) {
            return static_cast<Felem &>(*this) *= rhs.inv();
        }

        Felem operator/(const Felem &rhs) const {
            return Felem(static_cast<const Felem &>(*this)) /= rhs;
        }

        Integer getSize() const { return _f->getSize(); }

        const F getField() const { return *_f; }

        friend std::ostream &operator<<(std::ostream &os, const Felem &e) {
            os << to_string(e);
            return os;
        }

        friend bool operator==(Integer lhs, const Felem &rhs) {
            return rhs == lhs;
        }

        friend bool operator==(const Felem &lhs, Integer rhs) {
            return lhs == lhs.getField().get(rhs);
        }

        friend bool operator!=(Integer lhs, const Felem &rhs) {
            return !(lhs == rhs);
        }

        friend bool operator!=(const Felem &lhs, Integer rhs) {
            return !(lhs == rhs);
        }

        friend Felem &operator+=(Felem &lhs, Integer rhs) {
            lhs += lhs.getField().get(rhs);
            return lhs;
        }

        friend Felem operator+(const Felem &lhs, Integer rhs) {
            return lhs + lhs.getField().get(rhs);
        }

        friend Felem operator+(Integer lhs, const Felem &rhs) {
            return rhs.getField().get(lhs) + rhs;
        }

        friend Felem &operator-=(Felem &lhs, Integer rhs) {
            lhs -= lhs.getField().get(rhs);
            return lhs;
        }

        friend Felem operator-(const Felem &lhs, Integer rhs) {
            return lhs - lhs.getField().get(rhs);
        }

        friend Felem operator-(Integer lhs, const Felem &rhs) {
            return rhs.getField().get(lhs) - rhs;
        }

        friend Felem &operator*=(Felem &lhs, Integer rhs) {
            lhs *= lhs.getField().get(rhs);
            return lhs;
        }

        friend Felem operator*(const Felem &lhs, Integer rhs) {
            return lhs * lhs.getField().get(rhs);
        }

        friend Felem operator*(Integer lhs, const Felem &rhs) {
            return rhs.getField().get(lhs) * rhs;
        }

        friend bool compatible(const Felem &lhs, const Felem &rhs) {
            return lhs.getField() == rhs.getField();
        }

        friend Felem getZero(const Felem &e) { return e.getField().get(0); }

        friend Felem getOne(const Felem &e) { return e.getField().get(1); }

        friend std::string to_string(const Felem &e) {
            using std::to_string;
            return to_string(e._num);
        }

    protected:
        QuotientRing(Quotient num, F f) : _num(num), _mod(f.mod()), _f(new F(f)) { }

    private:
        bool initialized() const { return _f != nullptr; }

        void checkInSameField(const QuotientRing &rhs, std::string &&error) const {
            using std::to_string;
            if (this->getField() != rhs.getField())
                throw EOperationUnsupported(
                        error + "\nThe values that caused it were " + to_string(_num) +
                        " in " + to_string(this->getField()) +
                        " and " + to_string(rhs._num) +
                        " in F" + to_string(rhs.getField()));
        }

        Quotient _num;
        Quotient _mod;
        std::unique_ptr<F> _f;
    };
}

#endif // __FQELEM_HPP
