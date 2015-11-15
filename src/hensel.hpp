#ifndef __HENSEL
#define __HENSEL

#include<vector>
#include<utility>

class Zxelem;
class Fpxelem;

bool HenselLifting (const Zxelem &polynomial, unsigned int p, Fpxelem u1, Fpxelem w1, Zxelem & u, Zxelem & w);

#endif // __HENSEL
