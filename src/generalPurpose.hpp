#ifndef __GENERAL_PURPOSE_H
#define __GENERAL_PURPOSE_H

#include "types.hpp"

template<typename T, typename U> T fastPowMod(const T&a, U b, const T& p);
bool millerRabin(big_int n, int k = 35);
template<typename T> T eea (T a, T b, T &x, T &y);
template<typename T> T gcd(T a, T b);
template<typename T, typename U> T fastPow (const T& a, U b);

long long pollardRhoBrent (long long n);
bool pollardRhoLogarithm(long long g, long long h, long long n, long long & log);
#endif // __GENERAL_PURPOSE_H
