#ifndef __FELEMBASE_HPP
#define __FELEMBASE_HPP

#include "generalPurpose.hpp" // ExtendedEuclideanAlgorithm (eea)
#include "exceptions.hpp"

#include <stdexcept>        // std::runtime_error
#include <iosfwd>           // std::ostream
#include <memory>           // std::unique_ptr
#include <string>           // std::to_string


template<typename Felem, typename Fbase, typename Quotient, typename Integer>
class FelemBase{
    public:
        // Base field
        using F = Fbase ;

        FelemBase () : _f(nullptr){}
        FelemBase(const FelemBase & other) :
                _num(other._num),
                _mod(other._mod),
                _f (new F(*other._f))
                {}


        explicit operator Quotient() const { return _num; }

        Felem & operator=(const Felem &rhs){
            if(&rhs != this){
                if(!this->initialized()){
                    _num = rhs._num;
                    _mod = rhs._mod;
                    _f = std::unique_ptr<F>(new F(*rhs._f));
                }
                else{
                    checkInSameField(rhs, "Error in assignment.");
                    _num = rhs._num;
                }
            }
            return static_cast<Felem&>(*this);
        }

        Felem & operator=(Integer rhs){
            if(!this->initialized())
                throw std::runtime_error("Assignment to a non initialized Felem");
            static_cast<Felem&>(*this) = _f->get(rhs);
            return static_cast<Felem&>(*this);
        }

        bool operator==(const Felem &rhs)const{
            return (_num == rhs._num && *_f == *(rhs._f));
        }

        bool operator!=(const Felem &rhs)const{
            return !(static_cast<const Felem&>(*this) == rhs);
        }

        Felem & operator+=(const Felem &rhs){
            checkInSameField(rhs, "Addition or substraction error.");
            _num = (_num + rhs._num) % _mod;
            return static_cast<Felem&>(*this);
        }

        Felem operator+(const Felem &rhs) const{
            return Felem(static_cast<const Felem&>(*this)) += rhs;
        }

        Felem operator-() const{
            return _f->get(-_num);
        }

        Felem & operator-=(const Felem &rhs){
            return (static_cast<Felem&>(*this) +=(-rhs));
        }

        Felem operator-(const Felem &rhs) const{
            return Felem(static_cast<const Felem&>(*this)) -= rhs;
        }

        Felem & operator*=(const Felem &rhs){
            checkInSameField(rhs, "Multiplication or division error.");
            _num = (_num * rhs._num) % _mod;
            return static_cast<Felem&>(*this);
        }

        Felem operator*(const Felem &rhs) const{
            return Felem(static_cast<const Felem&>(*this)) *= rhs;
        }

        /** Multiplicative inverse */
        Felem inv() const{
            if(_num == 0)
                throw EOperationUnsupported("Error. Zero has no inverse.");
            Quotient res, aux;
            eea(_num, _mod, res, aux);
            return _f->get(res);
        }

        Felem & operator/=(const Felem &rhs){
            return static_cast<Felem&>(*this) *= rhs.inv();
        }

        Felem operator/(const Felem &rhs) const{
            return Felem(static_cast<const Felem&>(*this)) /= rhs;
        }

        Integer getSize() const{ return _f->getSize(); }

        const F getField() const{ return *_f; }

        friend std::ostream& operator<<(std::ostream& os, const Felem &e){
            os << to_string(e);
            return os;
        }

        friend bool operator==(Integer lhs, const Felem &rhs){
            return rhs == lhs;
        }
        friend bool operator==(const Felem &lhs, Integer rhs){
            return lhs == lhs.getField().get(rhs);
        }
        friend bool operator!=(Integer lhs, const Felem &rhs){
            return !(lhs == rhs);
        }
        friend bool operator!=(const Felem &lhs, Integer rhs){
            return !(lhs == rhs);
        }
        friend Felem & operator+=(Felem &lhs, Integer rhs){
            lhs+=lhs.getField().get(rhs);
            return lhs;
        }
        friend Felem operator+(const Felem &lhs, Integer rhs){
            return lhs + lhs.getField().get(rhs);
        }
        friend Felem operator+(Integer lhs, const Felem & rhs){
            return rhs.getField().get(lhs) + rhs;
        }
        friend Felem & operator-=(Felem &lhs, Integer rhs){
            lhs-=lhs.getField().get(rhs);
            return lhs;
        }
        friend Felem operator-(const Felem &lhs, Integer rhs){
            return lhs - lhs.getField().get(rhs);
        }
        friend Felem operator-(Integer lhs, const Felem & rhs){
            return rhs.getField().get(lhs) - rhs;
        }
        friend Felem & operator*=(Felem &lhs, Integer rhs){
            lhs*=lhs.getField().get(rhs);
            return lhs;
        }
        friend Felem operator*(const Felem &lhs, Integer rhs){
            return lhs * lhs.getField().get(rhs);
        }
        friend Felem operator*(Integer lhs, const Felem & rhs){
            return rhs.getField().get(lhs) * rhs;
        }

        friend bool compatible(const Felem &lhs, const Felem &rhs){
            return lhs.getField()==rhs.getField();
        }
        friend Felem getZero(const Felem &e){ return e.getField().get(0); }
        friend Felem getOne(const Felem &e){ return e.getField().get(1); }
        friend std::string to_string(const Felem &e){
            using std::to_string;
            return to_string(e._num);
        }

    protected:
        FelemBase (Quotient num, F f) : _num(num), _mod(f.mod()), _f(new F(f)){ }

    private:
        bool initialized() const{ return _f != nullptr; }
        void checkInSameField(const Felem &rhs, std::string&& error) const{
            using std::to_string;
            if(this->getField() != rhs.getField())
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

#endif // __FQELEM_HPP
