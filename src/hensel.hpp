#ifndef __HENSEL
#define __HENSEL

#include<vector>
#include<utility>
#include "henselSubsets.hpp"
#include "fpxelem.hpp"
#include "zxelem.hpp"


std::vector< std::pair < Zxelem_b, unsigned int > > squareFreeFactChar0(const Zxelem_b & pol);

bool HenselLifting (const Zxelem_b &polynomial, Fpxelem_b u1, Fpxelem_b w1, Zxelem_b & u, Zxelem_b & w);

std::vector< Zxelem_b > factorizationHenselSquareFree(Zxelem_b poli, HenselSubsets & hs);
std::vector< Zxelem_b > factorizationHenselSquareFree(const Zxelem_b & poli);
std::vector< std::pair < Zxelem_b, unsigned int > > factorizationHensel(const Zxelem_b & pol);
#endif // __HENSEL
