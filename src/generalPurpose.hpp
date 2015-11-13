#ifndef __GENERAL_PURPOSE_H
#define __GENERAL_PURPOSE_H

#include "types.hpp"

template<typename T, typename U> T fastPowMod(T a, U b, T p);
bool millerRabin(ll n, int k = 35);
template<typename T> T eea (T a, T b, T &x, T &y);
template<typename T> T gcd(T a, T b);
template<typename T> T fastPow (T a, ll b);
#endif // __GENERAL_PURPOSE_H
