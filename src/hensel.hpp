#ifndef __HENSEL
#define __HENSEL

#include<vector>
#include<utility>
#include "henselSubsets.hpp"

class Zxelem;
class Fpxelem;


bool HenselLifting (const Zxelem &polynomial, Fpxelem u1, Fpxelem w1, Zxelem & u, Zxelem & w);

std::vector< std::pair < Zxelem, unsigned int > > factorizationHensel(const Zxelem & pol);
std::vector< Zxelem > factorizationHenselSquareFree(const Zxelem & poli);
std::vector< Zxelem > factorizationHenselSquareFree(const Zxelem & poli, HenselSubsets & hs);
#endif // __HENSEL
