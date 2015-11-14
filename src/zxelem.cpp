#include "zxelem.hpp"
#include "fpxelem.hpp"

#include <vector>
#include <climits>

Zxelem::Zxelem(const big_int & e) : PolynomialRing<Zxelem, big_int>(e){}
Zxelem::Zxelem(const std::vector<big_int> & v) : PolynomialRing<Zxelem, big_int>(v){}
Zxelem::Zxelem(const Fpxelem & e) : PolynomialRing<Zxelem, big_int>(std::vector<big_int>(e._v.begin(), e._v.end())){}

Zxelem getZero(const Zxelem &e){return Zxelem(0);}
Zxelem getOne(const Zxelem &e){return Zxelem(1);}

const big_int unit(const Zxelem &e){
    return unit(e.lc());
}

bool compatible(const Zxelem &lhs, const Zxelem &rhs){ return true; }

big_int normInf(const Zxelem &e){
    big_int ret = e[0];
    for(int i=1;i<e.deg()+1;++i)
        if(e[i] > ret)
            ret = e[i];
    return ret;
}
