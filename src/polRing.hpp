#ifndef __POL_RING_HPP
#define __POL_RING_HPP

#include <vector>
#include <algorithm>        // find_if
#include <utility>          // pair, make_pair
#include <functional>       // not_equal_to
#include <string>           // to_string

#include "types.hpp"
#include "exceptions.hpp"

namespace alcp {
    template<template <class> class FxelemBase , class Felem, class Integer>
    class PolynomialRing {
    static_assert(is_integral<Integer>::value, "Type is not a supported integer.");
    private:
        using Fxelem = FxelemBase<Integer>;
    public:
        using Int = Integer;
        const static char var = 'x';

        PolynomialRing() = default;

        // Immersion from the base ring
        template<class Felem_t,
                class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>>
        PolynomialRing(const Felem_t &e) : _v(std::vector<Felem>({e})) { }

        // Move immersion from base ring
        template<class Felem_t,
                class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>>
        PolynomialRing(Felem_t &&e) : _v(std::vector<Felem>({std::move(e)})) { }

        PolynomialRing(const std::vector<Felem> &v) : _v(v) {
            if (v.size() != 0) {
                // Remove trailing zeros
                this->removeTrailingZeros();
                Felem aux = this->lc();
#ifndef ALCP_NO_CHECKS
                for (auto &e : v)
                    if (!compatible(aux, e))
                        throw ENotCompatible("Not all the elements in the array are in the same ring.");
#endif
            }
        }

        template<class Felem_t, class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>>
        PolynomialRing(const std::vector<Felem_t> &v) : _v(v) {
            if (v.size() != 0) {
                // Remove trailing zeros
                this->removeTrailingZeros();
                Felem aux = this->lc();
                for (auto &e : v)
                    if (!compatible(aux, e))
                        throw ENotCompatible("Not all the elements in the array are in the same ring.");
            }
        }

        PolynomialRing(const PolynomialRing &) = default;

        template<class Felem_t, class Int,
                class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>,
                class = std::enable_if_t<is_integral<Int>::value>>
        PolynomialRing(const PolynomialRing<FxelemBase, Felem_t, Int> & e)
            : _v(e._v.begin(), e._v.end()){};

        PolynomialRing(PolynomialRing &&) = default;

        template<class Felem_t, class Int,
                class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>,
                class = std::enable_if_t<is_integral<Int>::value>>
        PolynomialRing(PolynomialRing<FxelemBase, Felem_t, Int> && e)
            : _v(
                std::make_move_iterator(e._v.begin()),
                std::make_move_iterator(e._v.end())
                ){};

        // Copy assignment
        PolynomialRing &operator=(const PolynomialRing &rhs) {
            if (&rhs != this
#ifndef ALCP_NO_CHECKS
                    && rhs.init()) {
                if (this->init())
                    checkInSameField(rhs, "Assignation failed. The elements are not in the same ring.");
#else
                ){
#endif
                _v = rhs._v;
            }
            return *this;
        }

        template<class Felem_t, class Int,
                class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>,
                class = std::enable_if_t<is_integral<Int>::value>>
        PolynomialRing &operator=(const PolynomialRing<FxelemBase, Felem_t, Int> &rhs) {
            return *this = PolynomialRing(rhs);
        }

        // Move assignment
        PolynomialRing &operator=(PolynomialRing &&rhs) {
            if (&rhs != this
#ifndef ALCP_NO_CHECKS
                    && rhs.init()) {
                if (this->init())
                    checkInSameField(rhs, "Assignation failed. The elements are not in the same ring.");
#else
                ){
#endif
                _v = std::move(rhs._v);
            }
            return *this;
        }

        // Generalized move assignment
        template<class Felem_t, class Int,
                class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>,
                class = std::enable_if_t<is_integral<Int>::value>>
        PolynomialRing &operator=(PolynomialRing<FxelemBase, Felem_t, Int> &&rhs) {
            return *this = PolynomialRing(std::move(rhs));
        }

        // Copy assignment from base ring
        template<class Felem_t,
                class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>>
        Fxelem &operator=(const Felem_t &rhs) {
#ifndef ALCP_NO_CHECKS
            if (rhs.init() && this->init())
                checkInSameField(PolynomialRing(rhs),
                            "Assignation failed. The elements are not in the same ring.");
#endif

            _v = std::vector<Felem>(rhs);
            return static_cast<Fxelem&>(*this);
        }

        // Move assignment from base ring
        template<class Felem_t,
                class = std::enable_if_t<std::is_constructible<Felem, Felem_t>::value>>
        Fxelem &operator=(Felem_t &&rhs) {
#ifndef ALCP_NO_CHECKS
            if (rhs.init() && this->init())
                checkInSameField(PolynomialRing(rhs),
                            "Assignation failed. The elements are not in the same ring.");
#endif
            _v = std::vector<Felem>(std::move(rhs));
            return static_cast<Fxelem&>(*this);
        }

        explicit operator std::vector<Felem>() const { return _v; }

        friend inline bool operator==(const Fxelem &lhs, const Fxelem &rhs) {
            return lhs._v == rhs._v;
        }

        friend inline bool operator!=(const Fxelem &lhs, const Fxelem &rhs) {
            return lhs._v != rhs._v;
        }

        Fxelem &operator+=(const Fxelem &rhs) {
#ifndef ALCP_NO_CHECKS
            checkInSameField(PolynomialRing(rhs),
                        "Polynomials not in the same ring. Error when adding the polynomials.");
#endif

            auto v1 = _v.begin();
            auto v2 = rhs._v.begin();
            while (v1 != _v.end() && v2 != rhs._v.end()) {
                *v1 += *v2;
                ++v1;
                ++v2;
            }
            while (v2 != rhs._v.end()) {
                _v.push_back(*v2);
                ++v2;
            }
            this->removeTrailingZeros();
            return static_cast<Fxelem &>(*this);
        }

        friend inline Fxelem operator+(const Fxelem &lhs, const Fxelem &rhs) {
            return Fxelem(lhs) += rhs;
        }

        Fxelem operator-() const {
            std::vector<Felem> ret(_v.size(), getZero(this->lc()));

            for (size_t i = 0; i < _v.size(); ++i)
                ret[i] = Felem(-_v[i]);
            return Fxelem(ret);
        }

        Fxelem &operator-=(const Fxelem &rhs) {
            return (static_cast<Fxelem &>(*this) += (-rhs));
        }

        Fxelem operator-(const Fxelem &rhs) const{
            return Fxelem(static_cast<const Fxelem &>(*this)) -= rhs;
        }

        Fxelem &operator*=(const Fxelem &rhs) {
#ifndef ALCP_NO_CHECKS
            checkInSameField(PolynomialRing(rhs),
                        "Polynomials not in the same ring. Error when multiplying the polynomials.");
#endif
            std::vector<Felem> ret(rhs._v.size() + _v.size() - 1, getZero(this->lc()));
            for (size_t i = 0; i < _v.size(); ++i)
                for (size_t j = 0; j < rhs._v.size(); ++j) {
                    ret[i + j] += _v[i] * rhs._v[j];
                }
            _v = ret;
            this->removeTrailingZeros();
            return static_cast<Fxelem &>(*this);
        }

        friend inline Fxelem operator*(const Fxelem &lhs, const Fxelem &rhs) {
            return Fxelem(lhs) *= rhs;
        }

        // Implements long polynomial division
        // Return quotient and reminder in first and second respectively
        std::pair<Fxelem, Fxelem> div2(const Fxelem &divisor) const {
#ifndef ALCP_NO_CHECKS
            checkInSameField(PolynomialRing(divisor),
                        "Polynomials not in the same ring. Error when dividing the polynomials.");
#endif

            if (divisor == 0)
                throw EOperationUnsupported("Error. Cannot divide by the polynomial 0");
            if (this->deg() < divisor.deg())
                return std::make_pair(Fxelem(getZero(this->lc())), static_cast<const Fxelem &>(*this));

            if (divisor.deg() == 0) {
                Fxelem quot(static_cast<const Fxelem &>(*this));
                Fxelem rem(static_cast<const Fxelem &>(*this));
                for (auto &e : quot._v) {
                    e /= divisor.lc();
                }
                return std::make_pair(quot, Fxelem(getZero(this->lc())));
            }

            Fxelem quot(getZero(this->lc()));
            Fxelem rem(static_cast<const Fxelem &>(*this));

            while (rem.deg() >= divisor.deg()) {
                std::vector<Felem> paddingZeros(rem.deg() - divisor.deg(), getZero(this->lc()));

                paddingZeros.push_back(rem.lc() / divisor.lc());
                Fxelem monDiv(paddingZeros);
                quot += monDiv;
                rem -= monDiv * divisor;
            }
            rem.removeTrailingZeros();

            return std::make_pair(quot, rem);
        }

        Fxelem &operator/=(const Fxelem &rhs) {
            *this = this->div2(rhs).first;
            return static_cast<Fxelem &>(*this);
        }

        Fxelem operator/(const Fxelem &rhs) const {
            return Fxelem(static_cast<const Fxelem &>(*this)) /= rhs;
        }

        Fxelem &operator%=(const Fxelem &rhs) {
            *this = this->div2(rhs).second;
            return static_cast<Fxelem &>(*this);
        }

        Fxelem operator%(const Fxelem &rhs) const {
            return Fxelem(static_cast<const Fxelem &>(*this)) %= rhs;
        }

        const Felem &operator[](size_t i) const { return _v[i]; }

        Felem &operator[](size_t i) { return _v[i]; }

        Fxelem derivative() const {
            if (this->deg() == 0)
                return Fxelem(getZero(this->lc()));
            std::vector<Felem> v(_v);
            for (size_t i = 1; i < v.size(); ++i)
                v[i - 1] = v[i] * i;
            v.pop_back();
            Fxelem ret(v);
            ret.removeTrailingZeros();

            return ret;
        }

        // Leading coefficient
        Felem lc() const { return _v.back(); }

        // Degree of the polynomial
        size_t deg() const { return _v.size() - 1; }

        std::size_t nonZeroCoefs() const{
            std::size_t ret = 0;
            for(const auto& elem : _v)
                if(elem != 0)
                    ++ret;
            return ret;
        }

        // Normal form of the polynomial. It ensures the unicity of gdc for example
        friend inline Fxelem normalForm(const Fxelem &e) { return e / unit(e); }

        friend std::string to_string(const Fxelem &f, char var = PolynomialRing::var /* x */) {
            if(f.deg() == 0)
                return to_string(f.lc());

            std::string aux(to_string_coef(f.lc()));
            if(aux[0] == '+'){
                aux.erase(0,1);
                if(aux == "1")
                    aux = "";
            }
            else if(aux == "-1")
                aux = "-";
            std::string s(aux);

            s += var;
            if(f.deg() != 1)
                s += "^" + std::to_string(f.deg());

            for (int i = static_cast<int>(f._v.size()) - 2; i >= 1; --i) {
                if (f._v[i] == 0)
                    continue;
                aux = to_string_coef(f._v[i]);
                if (aux == "+1")
                    s += "+";
                else if (aux == "-1")
                    s += "-";
                else
                    s += aux;

                s += var;
                if(i != 1)
                    s += "^" + std::to_string(i);

            }
            if(f._v[0] != 0)
                s += to_string_coef(f._v[0]);
            return s;
        }

        friend inline std::ostream& operator<<(std::ostream &os, const Fxelem &f) {
            return os << to_string(f);
        }

    protected:
        std::vector<Felem> _v;

    private:

        template <template <class> class, class, class>
        friend class PolynomialRing;

        void removeTrailingZeros() {
            Felem zero = getZero(this->lc());
            _v.erase(
                    std::find_if(
                            _v.rbegin(),
                            _v.rend(),
                            std::bind1st(std::not_equal_to<Felem>(), zero)).base(),
                    _v.end());
            // In case it was the polynomial equal to zero
            if (_v.size() == 0)
                _v.push_back(std::move(zero));
        }

#ifndef ALCP_NO_CHECKS
        // Relies in the fact that it is not possible to quotient by the ideal generated by 0
        bool init() const { return _v.size() != 0; }

        void checkInSameField(const PolynomialRing &rhs, std::string &&error) const {
            if (!compatible(this->lc(), rhs.lc()))
                throw EOperationUnsupported(
                        error +
                        "\nThe values that caused it were " +
                        to_string(static_cast<const Fxelem &>(*this)) +
                        " and " +
                        to_string(static_cast<const Fxelem &>(rhs)) +
                        " are not in the same ring.");
        }
#endif

    };
}

#endif //  __POL_RING_HPP
