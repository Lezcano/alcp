#ifndef __FELEMBASE_HPP
#define __FELEMBASE_HPP

#include <iosfwd>           // std::ostream
#include <memory>           // std::unique_ptr,  std::make_unique
#include <string>           // std::to_string
#include <type_traits>      // std::enable_if, std::is_integral
#include <utility>          // std::move

#include "generalPurpose.hpp" // ExtendedEuclideanAlgorithm (eea)
#include "exceptions.hpp"

namespace alcp {
    template<template <class> class FelemBase, class Quotient, class Integer>
    class QuotientRing {
    static_assert(is_integral<Integer>::value, "Type is not a supported integer.");
    public:
        using Felem = FelemBase<Integer>;
        using Int = Integer;
        QuotientRing() = default;

        QuotientRing(const QuotientRing &) = default;

        // Necesario Integer porque si hacemos class Felem, en vez de class template,
        //  no podemos hacer const QuotientRing<Felem,Quot> con Quot un template dado que
        //  los templates serian diferentes
        template<class Int, class Quot,
                class = std::enable_if_t<is_integral<Int>::value>>
        QuotientRing(const QuotientRing<FelemBase, Quot, Int> &rhs) :_mod(rhs._mod), _num(rhs._num){ }

        QuotientRing(QuotientRing &&) = default;

        template<class Int, class Quot,
                class = std::enable_if_t<is_integral<Int>::value>>
        QuotientRing(QuotientRing<FelemBase, Quot, Int> &&rhs) : _mod(std::move(rhs._mod)), _num(std::move(rhs._num)){ }

        // Copy assignment
#ifndef ALCP_NO_CHECKS
        QuotientRing &operator=(const QuotientRing &rhs) {
            if (&rhs != this && rhs.init()){
                if (this->init())
                    checkInSameField(rhs, "The elements are not in the same ring.");
                _mod = rhs._mod;
                _num = rhs._num;
            }
            return *this;
        }
#else
        QuotientRing &operator=(const QuotientRing &rhs) = default;
#endif

        // Duplicated code for generalized copy... :/
        template<class Int, class Quot,
                class = std::enable_if_t<is_integral<Int>::value>>
        QuotientRing &operator=(const QuotientRing<FelemBase, Quot, Int> &rhs){
            return *this = QuotientRing(rhs);
        }

        // Move assignment
#ifndef ALCP_NO_CHECKS
        QuotientRing &operator=(QuotientRing &&rhs) {
            if (&rhs != this && rhs.init()) {
                if (this->init())
                    checkInSameField(rhs, "The elements are not in the same ring.");
                _mod = std::move(rhs._mod);
                _num = std::move(rhs._num);
            }
            return *this;
        }
#else
        QuotientRing &operator=(QuotientRing &&rhs) = default;
#endif

        template<class Int, class Quot,
                class = std::enable_if_t<is_integral<Int>::value>>
        QuotientRing &operator=(QuotientRing<FelemBase, Quot, Int> &&rhs){
            return *this = QuotientRing(std::move(rhs));
        }

        template<class Int2, class = typename std::enable_if<is_integral<Int2>::value>::type>
        Felem &operator=(Int2 rhs) {
#ifndef ALCP_NO_CHECKS
            if (!this->init())
                throw ENotInitializedRing("Assignment to a non initialized QuotientRing elem.");
#endif
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

        friend inline bool operator<(const Felem &lhs, const Felem &rhs){
			return (lhs._num < rhs._num && lhs._mod == rhs._mod);
		}

        friend inline bool operator<=(const Felem &lhs, const Felem &rhs){
			return (lhs._num <= rhs._num && lhs._mod == rhs._mod);
		}

        friend inline bool operator>(const Felem &lhs, const Felem &rhs){
			return (lhs._num > rhs._num && lhs._mod == rhs._mod);
		}

		friend inline bool operator>=(const Felem &lhs, const Felem &rhs){
			return (lhs._num >= rhs._num && lhs._mod == rhs._mod);
		}

        Felem &operator+=(const Felem &rhs) {
#ifndef ALCP_NO_CHECKS
            checkInSameField(rhs, "Addition or substraction error.");
#endif
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
#ifndef ALCP_NO_CHECKS
            checkInSameField(rhs, "Multiplication or division error.");
#endif
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

        auto getSize() const { return static_cast<const Felem*>(this)->getField().getSize(); }


        friend std::ostream &operator<<(std::ostream &os, const Felem &e) {
            os << to_string(e);
            return os;
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend bool operator==(const Felem &lhs, Int rhs) {
            return lhs == lhs.getField().get(rhs);
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend bool operator==(Int lhs, const Felem &rhs) {
            return rhs == lhs;
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend bool operator!=(const Felem &lhs, Int rhs) {
            return !(lhs == rhs);
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend bool operator!=(Int lhs, const Felem &rhs) {
            return !(lhs == rhs);
        }


        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem &operator+=(Felem &lhs, Int rhs) {
            return lhs += lhs.getField().get(rhs);
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem operator+(const Felem &lhs, Int rhs) {
            return lhs + lhs.getField().get(rhs);
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem operator+(Int lhs, const Felem &rhs) {
            return rhs.getField().get(lhs) + rhs;
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem &operator-=(Felem &lhs, Int rhs) {
            return lhs -= lhs.getField().get(rhs);
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem operator-(const Felem &lhs, Int rhs) {
            return lhs - lhs.getField().get(rhs);
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem operator-(Int lhs, const Felem &rhs) {
            return rhs.getField().get(lhs) - rhs;
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem &operator*=(Felem &lhs, Int rhs) {
            return lhs *= lhs.getField().get(rhs);
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem operator*(const Felem &lhs, Int rhs) {
            return lhs * lhs.getField().get(rhs);
        }

        template<class Int, class = std::enable_if_t<is_integral<Int>::value>>
        friend Felem operator*(Int lhs, const Felem &rhs) {
            return rhs.getField().get(lhs) * rhs;
        }

        friend bool compatible(const Felem &lhs, const Felem &rhs) {
            return lhs.getField() == rhs.getField();
        }

        friend Felem getZero(const Felem &e) { return e.getField().get(0); }

        friend Felem getOne(const Felem &e) { return e.getField().get(1); }

    protected:
        QuotientRing(const Quotient& num, const Quotient& mod) : _num(num), _mod(mod){ }
        QuotientRing(Quotient&& num, Quotient&& mod) : _num(std::move(num)), _mod(std::move(mod)){ }

        Quotient _num;
        Quotient _mod;

        template <template <class> class, class, class>
            friend class QuotientRing;

        // It shows the init function
        template <template <class> class, class, class>
            friend class PolynomialRing;
    private:

#ifndef ALCP_NO_CHECKS
        // Relies in the fact that it is not possible to quotient by the ideal generated by 0
        bool init() const { return _mod != Quotient(); }

        void checkInSameField(const QuotientRing &rhs, std::string &&error) const {
            if (!compatible(static_cast<const Felem &>(*this),
                            static_cast<const Felem &>(rhs)))
                throw EOperationUnsupported(
                        error +
                        "\nThe values that caused it were " +
                        to_string(_num) +
                        " in " +
                        to_string(static_cast<const Felem*>(this)->getField()) +
                        " and " +
                        to_string(rhs._num) +
                        " in " +
                        to_string(static_cast<const Felem &>(rhs).getField()));
        }
#endif
    };
}

#endif // __FELEMBASE_HPP
