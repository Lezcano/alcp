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
            class Quotient, class Integer>
    class QuotientRing {
    public:
        // Base field
        using Felem = FelemBase<Integer>;

        QuotientRing() = default;

        QuotientRing(const QuotientRing &other) = default;

        QuotientRing(QuotientRing &&) = default;

        // Generalized copy constructor
        template<class Int2, class = typename std::enable_if<std::is_integral<Int2>::value>::type>
        QuotientRing(const QuotientRing<FelemBase, Quotient, Int2> &other) :
                _num(other._num),
                _mod(other._mod),
                _init(other._init)
        { }

        // Copy assignment
        QuotientRing &operator=(const QuotientRing &rhs) {
            if (&rhs != this && rhs._init) {
                if (this->_init)
                    checkInSameField(rhs, "The elements are not in the same ring.");
                _init = true;
                _mod = rhs._mod;
                _num = rhs._num;
            }
            return static_cast<Felem &>(*this);
        }

        // Generalized copy assignment
        template<class Int2, class = typename std::enable_if<std::is_integral<Int2>::value>::type>
        QuotientRing &operator=(const QuotientRing<FelemBase, Quotient, Int2> &rhs) {
            if (&rhs != this && rhs._init) {
                if (this->_init)
                    checkInSameField(rhs, "The elements are not in the same ring.");
                _init = true;
                _mod = rhs._mod;
                _num = rhs._num;
            }
            return static_cast<Felem &>(*this);
        }

        // Generalized move assignment
        template<class Int2, class = typename std::enable_if<std::is_integral<Int2>::value>::type>
        QuotientRing &operator=(QuotientRing<FelemBase, Quotient, Int2> &&rhs) {
            if (&rhs != this && rhs._init) {
                if (this->_init)
                    checkInSameField(rhs, "The elements are not in the same ring.");
                _init = true;
                _mod = std::move(rhs._mod);
                _num = std::move(rhs._num);
            }
            return static_cast<Felem &>(*this);
        }

        template<class Int2, class = typename std::enable_if<std::is_integral<Int2>::value>::type>
        Felem &operator=(Int2 rhs) {
            if (!this->_init)
                throw std::runtime_error("Assignment to a non initialized QuotientRing elem.");
            static_cast<Felem &>(*this) =
                static_cast<const Felem*>(this)->getField().get(rhs);
            return static_cast<Felem &>(*this);
        }

        explicit operator Quotient() const { return _num; }

        friend inline bool operator==(const Felem &lhs, const Felem &rhs){
            return (lhs._num == rhs._num && lhs._mod == rhs._mod);
        }

        friend inline bool operator!=(const Felem &lhs, const Felem &rhs){
            return !(lhs == rhs);
        }

        Felem &operator+=(const Felem &rhs) {
            checkInSameField(rhs, "Addition or substraction error.");
            _num = (_num + rhs._num) % _mod;
            return static_cast<Felem &>(*this);
        }

        friend inline Felem operator+(const Felem& lhs, const Felem &rhs){
            return Felem(lhs) += rhs;
        }

        Felem operator-() const {
            return static_cast<const Felem*>(this)->getField().get(-_num);
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

        friend inline Felem operator*(const Felem &lhs, const Felem &rhs){
            return Felem(lhs) *= rhs;
        }

        /** Multiplicative inverse */
        Felem inv() const {
            if (_num == 0)
                throw EOperationUnsupported("Error. Zero has no inverse.");
            Quotient res, aux;
            eea(_num, _mod, res, aux);

            return static_cast<const Felem*>(this)->getField().get(res);
        }

        Felem &operator/=(const Felem &rhs) {
            return static_cast<Felem &>(*this) *= rhs.inv();
        }

        Felem operator/(const Felem &rhs) const {
            return Felem(static_cast<const Felem &>(*this)) /= rhs;
        }

        Integer getSize() const { return static_cast<Felem*>(this)->getField().getSize(); }


        friend std::ostream &operator<<(std::ostream &os, const Felem &e) {
            os << to_string(e);
            return os;
        }

        friend bool operator==(const Felem &lhs, Integer rhs) {
            return lhs == lhs.getField().get(rhs);
        }

        friend bool operator==(Integer lhs, const Felem &rhs) {
            return rhs == lhs;
        }

        friend bool operator!=(const Felem &lhs, Integer rhs) {
            return !(lhs == rhs);
        }

        friend bool operator!=(Integer lhs, const Felem &rhs) {
            return !(lhs == rhs);
        }


        friend Felem &operator+=(Felem &lhs, Integer rhs) {
            return lhs += lhs.getField().get(rhs);
        }

        friend Felem operator+(const Felem &lhs, Integer rhs) {
            return lhs + lhs.getField().get(rhs);
        }

        friend Felem operator+(Integer lhs, const Felem &rhs) {
            return rhs.getField().get(lhs) + rhs;
        }

        friend Felem &operator-=(Felem &lhs, Integer rhs) {
            return lhs -= lhs.getField().get(rhs);
        }

        friend Felem operator-(const Felem &lhs, Integer rhs) {
            return lhs - lhs.getField().get(rhs);
        }

        friend Felem operator-(Integer lhs, const Felem &rhs) {
            return rhs.getField().get(lhs) - rhs;
        }

        friend Felem &operator*=(Felem &lhs, Integer rhs) {
            return lhs *= lhs.getField().get(rhs);
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
        QuotientRing(const Quotient& num, const Quotient& mod) : _num(num), _mod(mod), _init(true){ }

        Quotient _num;
        Quotient _mod;

    private:

        void checkInSameField(const QuotientRing &rhs, std::string &&error) const {
            using std::to_string;
            if (!compatible(static_cast<const Felem &>(*this),
                            static_cast<const Felem &>(rhs)))
                throw EOperationUnsupported(
                error + "\nThe values that caused it were " +
                to_string(_num) + " in " +
                to_string(static_cast<const Felem*>(this)->getField()) +
                " and " + to_string(rhs._num) + " in F" +
                to_string(static_cast<const Felem &>(rhs).getField()));
        }

        bool _init{false};
    };
}

#endif // __FQELEM_HPP
