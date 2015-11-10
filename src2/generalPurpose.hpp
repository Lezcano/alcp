#ifndef __GENERAL_PURPOSE_H
#define __GENERAL_PURPOSE_H

#include "types.hpp"

template<typename T> T fastPowMod(T a, ll b, T p);
bool millerRabin(ll n, int k = 35);
template<typename T> T eea (T a, T b, T &x, T &y);
template<typename T> T gcd(T a, T b);

#endif // __GENERAL_PURPOSE_H
