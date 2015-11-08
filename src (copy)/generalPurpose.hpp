#ifndef __GENERAL_PURPOSE_H
#define __GENERAL_PURPOSE_H

#include "types.hpp"

ll fastPowMod(ll a, ll b, ll p);
bool millerRabin(ll n, int k = 35);
ll eea (ll a, ll b, ll& x, ll& y);
ll gcd(ll a, ll b);
template<typename T> T eea (T a, T b, T &x, T &y);
template<typename T> T gcd(T a, T b);

#endif // __GENERAL_PURPOSE_H
