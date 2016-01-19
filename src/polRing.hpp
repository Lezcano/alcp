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
    template<typename Fxelem, typename Felem>
    class PolynomialRing {
    public:

        const static char var = 'x';

        PolynomialRing() = default;

        // Immersion from the base ring
        PolynomialRing(const Felem &e) : _v(std::vector<Felem>({e})) { }

        PolynomialRing(const std::vector<Felem> &v) : _v(v) {
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
        PolynomialRing(PolynomialRing &&) = default;

        // Copy assignment
        PolynomialRing &operator=(const PolynomialRing &rhs) {
            if (&rhs != this && rhs.initialized()) {
                if (this->initialized() && !compatible(this->lc(), rhs.lc()))
                    throw ENotCompatible(
                            "Asignation failed. The vectors " +
                                to_string(static_cast<Fxelem &>(*this)) +
                                " and " +
                                to_string(static_cast<const Fxelem &>(rhs)) +
                                " are not in the same ring.");
                _v = rhs._v;
            }
            return *this;
        }

        // Move assignment
        PolynomialRing &operator=(PolynomialRing &&rhs) {
            if (&rhs != this && rhs.initialized()) {
                if (this->initialized() && !compatible(this->lc(), rhs.lc()))
                    throw ENotCompatible(
                            "Asignation failed. The vectors " +
                                to_string(static_cast<Fxelem &>(*this)) +
                                " and " +
                                to_string(static_cast<Fxelem &&>(rhs)) +
                                " are not in the same ring.");
                _v = std::move(rhs._v);
            }
            return *this;
        }

        // Copy assignment from base ring
        Fxelem &operator=(const Felem &rhs) {
            if (this->initialized() && !compatible(this->lc(), rhs))
                throw ENotCompatible(
                        "Asignation failed. The vector " +
                            to_string(static_cast<Fxelem &>(*this)) +
                            " and the element " +
                            to_string(static_cast<const Fxelem &>(rhs)) +
                            " are not in the same ring.");
            _v = std::vector<Felem>(rhs);
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
            if (!compatible(static_cast<Fxelem &>(*this), rhs))
                throw EOperationUnsupported(
                        "Polynomials not in the same ring. Error when adding the polynomials " +
                        to_string(static_cast<Fxelem &>(*this)) + " and " + to_string(rhs) + ".");
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
            if (!compatible(static_cast<Fxelem &>(*this), rhs))
                throw EOperationUnsupported(
                        "Polynomials not in the same ring. Error when multiplying the polynomials " +
                        to_string(static_cast<Fxelem &>(*this)) +
                        " and " + to_string(rhs) + ".");

            std::vector<Felem> ret(rhs._v.size() + _v.size() - 1, getZero(this->lc()));
            for (size_t i = 0; i < _v.size(); ++i)
                for (size_t j = 0; j < rhs._v.size(); ++j) {
                    ret[i + j] += _v[i] * rhs._v[j];
                }
            _v = ret;
            return static_cast<Fxelem &>(*this);
        }

        friend inline Fxelem operator*(const Fxelem &lhs, const Fxelem &rhs) {
            return Fxelem(lhs) *= rhs;
        }

        // Implements long polynomial division
        // Return quotient and reminder in first and second respectively
        std::pair<Fxelem, Fxelem> div2(const Fxelem &divisor) const {
            if (!compatible(static_cast<const Fxelem &>(*this), divisor))
                throw EOperationUnsupported(
                        "Polynomials not in the same ring. Error when dividing the polynomials " +
                        to_string(static_cast<const Fxelem &>(*this)) +
                        " and " + to_string(divisor) + ".");

            if (divisor.deg() == 0 && divisor._v[0] == 0)
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


        friend inline Fxelem operator+(const Fxelem &lhs, const Felem &rhs){ return lhs + Fxelem(rhs); }
        friend inline Fxelem operator-(const Fxelem &lhs, const Felem &rhs){ return lhs - Fxelem(rhs); }
        friend inline Fxelem operator*(const Fxelem &lhs, const Felem &rhs){ return lhs * Fxelem(rhs); }
        friend inline Fxelem operator/(const Fxelem &lhs, const Felem &rhs){ return lhs / Fxelem(rhs); }
        friend inline Fxelem operator%(const Fxelem &lhs, const Felem &rhs){ return lhs % Fxelem(rhs); }
        friend inline bool operator==(const Fxelem &lhs, const Felem &rhs){ return lhs == Fxelem(rhs); }
        friend inline bool operator!=(const Fxelem &lhs, const Felem &rhs){ return lhs != Fxelem(rhs); }

        friend inline Fxelem &operator+=(Fxelem &lhs, const Felem &rhs){return lhs += Fxelem(rhs); }
        friend inline Fxelem &operator-=(Fxelem &lhs, const Felem &rhs){ return lhs -= Fxelem(rhs); }
        friend inline Fxelem &operator*=(Fxelem &lhs, const Felem &rhs){ return lhs *= Fxelem(rhs); }
        friend inline Fxelem &operator/=(Fxelem &lhs, const Felem &rhs){ return lhs /= Fxelem(rhs); }
        friend inline Fxelem &operator%=(Fxelem &lhs, const Felem &rhs){ return lhs %= Fxelem(rhs); }

        friend inline Fxelem operator+(const Felem &lhs, const Fxelem &rhs){ return Fxelem(lhs) + rhs; }
        friend inline Fxelem operator-(const Felem &lhs, const Fxelem &rhs){ return Fxelem(lhs) - rhs; }
        friend inline Fxelem operator*(const Felem &lhs, const Fxelem &rhs){ return Fxelem(lhs) * rhs; }
        friend inline Fxelem operator/(const Felem &lhs, const Fxelem &rhs){ return Fxelem(lhs) / rhs; }
        friend inline Fxelem operator%(const Felem &lhs, const Fxelem &rhs){ return Fxelem(lhs) % rhs; }
        friend inline bool operator==(const Felem &lhs, const Fxelem &rhs){ return Fxelem(lhs) == rhs; }
        friend inline bool operator!=(const Felem &lhs, const Fxelem &rhs){ return Fxelem(lhs) != rhs; }

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

        bool initialized() const { return _v.size() != 0; }

    };
}

#endif //  __POL_RING_HPP
