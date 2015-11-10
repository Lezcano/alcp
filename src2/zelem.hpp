#ifndef __ZELEM_HPP
#define __ZELEM_HPP

#include "types.hpp"
#include "polRing.hpp"

#include <string>

std::string to_string(bint e);
bool compatible(bint lhs, bint rhs);
const bint unit(bint e);
const bint normalForm(bint e);
bint getZero(bint e);
bint getOne(bint e);



class Zxelem : public PolinomialRing<Zxelem, bint>{
    public:
        Zxelem(const bint & e);
        Zxelem(const std::vector<bint> & v);
};

const bint unit(const Zxelem &e);

#endif
#endif // __ZELEM_HPP
