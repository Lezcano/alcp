#ifndef __HENSEL
#define __HENSEL

#include<vector>
#include<utility>
#include "henselSubsets.hpp"
#include "fpxelem.hpp"
#include "zxelem.hpp"


bool HenselLifting (const Zxelem_b &polynomial, Fpxelem_b u1, Fpxelem_b w1, Zxelem_b & u, Zxelem_b & w);

std::vector< std::pair < Zxelem_b, unsigned int > > factorizationHensel(const Zxelem_b & pol);
std::vector< Zxelem_b > factorizationHenselSquareFree(const Zxelem_b & poli);
std::vector< Zxelem_b > factorizationHenselSquareFree(const Zxelem_b & poli, HenselSubsets & hs);
#endif // __HENSEL
