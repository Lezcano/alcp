#include "zxelem.hpp"
#include "fpxelem.hpp"
#include "generalPurpose.hpp" // gcd

#include <vector>
#include <algorithm> // transform

// Auxiliary function for the Zxelem ctor
// It returns the symmetric representation of f \in Fp[X] as an element of Z[X]
Zxelem toZxelemSym(const Fpxelem &e){
    std::vector<big_int> v(e._v.size());
    big_int p = e.getSize();
    // When p = 2, p2 = 1 and so every v[i] is <= p2
    // When p is odd, p/2 = (p-1)/2 and so we get the symmetric representation
    big_int p2 = p/2;
    std::transform(e._v.begin(), e._v.end(), v.begin(),
            [&p, &p2](const Fpelem &e) -> big_int {
            return static_cast<big_int>(e) <= p2 ?
                static_cast<big_int>(e) :
                static_cast<big_int>(e)-static_cast<big_int>(p);
                });
    return Zxelem(v);
}

Zxelem::Zxelem() = default;
Zxelem::Zxelem(const big_int & e) : PolynomialRing<Zxelem, big_int>(e){}
Zxelem::Zxelem(const std::vector<big_int> & v) : PolynomialRing<Zxelem, big_int>(v){}
Zxelem::Zxelem(const Fpxelem & e) : PolynomialRing<Zxelem, big_int>(toZxelemSym(e)){}

Zxelem getZero(const Zxelem &e){return Zxelem(0);}
Zxelem getOne(const Zxelem &e){return Zxelem(1);}

big_int unit(const Zxelem &e){
    return unit(e.lc());
}

bool compatible(const Zxelem &lhs, const Zxelem &rhs){ return true; }

big_int normInf(const Zxelem &e){
    big_int ret = e[0];
    for(size_t i=1;i<=e.deg();++i)
        if(e[i] > ret)
            ret = e[i];
    return ret;
}

big_int content(const Zxelem & e){
    big_int gcdE = e[0];
    for(size_t i=1;i<=e.deg();++i)
        gcdE = gcd(gcdE, e[i]);
    return gcdE;
}
