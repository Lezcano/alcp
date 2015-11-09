#include <vector>
#include <string>

#include "fq.hpp"
// types.hpp included in fp.hpp
#include "fqelem.hpp"
#include "exceptions.hpp"
#include "generalPurpose.hpp" // Miller Rabin

Fq::Fq(ll p, int n): _p(p), _n(n), _base(p){
    std::vector<Fpelem> _mod (_n, _base.get(0));

    v.back() = _base.get(1);
    v[0] = _base.get(1);
    // Hopefully this ends in a reasonable amount of time
    //  Maybe set a timer?
    //  The probability of getting an irreducible polinomial
    //   in each iteration is 1/n
    while(!Fpxelem(v).irreducible())
        increment(v);

    _mod = Fpxelem(v);
}

Fqelem Fq::get(ll n)const{
    return Fqelem(Fpxelem(n), _mod, std::unique_ptr<Fq>(new Fq(*this)));
}

Fqelem Fq::get(Fpxelem f)const{
    return Fqelem(f, _mod, std::unique_ptr<Fq>(new Fq(*this)));
}

ll Fq::getSize()const{return (ll) std::pow((double)_p,(double)n);}

bool Fq::operator==(const Fq &rhs)const{return _p == rhs._p && _n == rhs._n;}
bool Fq::operator!=(const Fq &rhs)const{return _p != rhs._p && _n != rhs._n;}

std::vector<Fqelem> Fq::getElems(){
    std::vector<Fqelem> ret;
    std::vector<Fpelem> act(_n, _base.get(0));

    do {
        ret.push_back(this->get(Fpxelem(act)));
    } while (increase(act));
    return ret;
}

std::string to_string(const Fq &e){
    return "GF(" + std::to_string(e._p) + "^" std::to_string(e._n) + ")";
}

bool increase(std::vector<Fpelem>& act) {
    for(auto e = act.begin(); e != act.end(); ++e){
        *e += 1;
        if (*e != 0)
            return true;
    }
    return false; // We are done
}
