#ifndef __ZXELEM_HPP
#define __ZXELEM_HPP

#include "types.hpp"
#include "zelem.hpp"    // Auxliary functions
#include "fpxelem.hpp"
#include "polRing.hpp"

template<class Integer> class Fpxelem;

template<class Integer>
class Zxelem : public PolynomialRing<Zxelem<Integer>, Integer>{
    public:
        using PolynomialRing<Zxelem<Integer>, Integer>::PolynomialRing;
        Zxelem() = default;

		/* It returns a polynomial in Zx using the symmetric
		 * representation of Fp, this is
		 * -(p-1)/2, -(p-1)/2+1, ..-1, 0, 1, .. (p-1)/2
		 * */
        Zxelem(const Fpxelem<Integer> & e) : PolynomialRing<Zxelem<Integer>, Integer>(toZxelemSym(e)){}

        friend class Fpxelem<Integer>;
        friend Fpxelem<Integer> toFpxelem(const Zxelem<Integer> &e, Integer p);

        friend Zxelem<Integer> getZero(const Zxelem<Integer> &e){ return Zxelem<Integer>(0); }
        friend Zxelem<Integer> getOne(const Zxelem<Integer> &e){ return Zxelem<Integer>(1); }
        friend Integer unit(const Zxelem<Integer> &e){ return unit(e.lc()); }
        friend bool compatible(const Zxelem<Integer> &lhs, const Zxelem<Integer> &rhs){ return true; }
        friend Integer normInf(const Zxelem<Integer> &e){
            Integer ret = e[0];
            for(std::size_t i=1;i<=e.deg();++i)
                if(e[i] > ret)
                    ret = e[i];
            return ret;
        }
        friend Integer content(const Zxelem<Integer> &e){
            Integer gcdE = e[0];
            for(size_t i = 1;i <= e.deg();++i)
                gcdE = gcd(gcdE, e[i]);
            return gcdE;
        }

    private:

        // It returns the symmetric representation of f \in Fp[X] as an element of Z[X]
        Zxelem toZxelemSym(const Fpxelem<Integer> &e){
            std::vector<Integer> v(e._v.size());
            Integer p = e.getSize();
            // When p = 2, p2 = 1 and so every v[i] is <= p2
            // When p is odd, p/2 = (p-1)/2 and so we get the symmetric representation
            Integer p2 = p/2;
            std::transform(e._v.begin(), e._v.end(), v.begin(),
                    [&p, &p2](const Fpelem<Integer> &e) -> Integer {
                    return static_cast<Integer>(e) <= p2 ?
                        static_cast<Integer>(e) :
                        static_cast<Integer>(e)-static_cast<Integer>(p);
                        });
            return Zxelem(v);
        }
};

using Zxelem_b = Zxelem<big_int>;

#endif
