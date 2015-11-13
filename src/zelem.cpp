#include <string>
#include <sstream>

#include "zelem.hpp"
// Types inlcuded in zelem.hpp

/**
 * Addition of 63 bit numbers without overflow
 */
long long add(long long a, long long b, long long p){
    return (long long)(((unsigned long long) a + (unsigned long long) b) % p);
}

/** Russian peasant multiplication
 *   Overview
 *    It multiplies two positive integers of 63 bits and reduces them
 *     modulo p, using integers not bigger than 64 bits.
 *    It circumvents the problem of not having integers greater
 *     than 64 bits in C++.
 *    It does this by computing the multiplication adding 2^i*b(mod p)
 *     to the result if the i-th bit of a is one.
 */
long long russianPeasantMultiplication(long long a, long long b, long long p){
    unsigned long long res = 0;
    while (a != 0) {
        if (a & 1) res = (res + b) % p;
        a >>= 1;
        b = (b << 1) % p;
    }
    return (long long) res;
}
bool compatible(big_int lhs, big_int rhs){
    return true;
}
const big_int unit(big_int e){
    return e >= 0 ? (big_int)1 : (big_int)-1;
}
const big_int normalForm(big_int e){ return e/unit(e);}

//std::string to_string(big_int e){return e.convert_to<std::string>();}
std::string to_string(big_int e){return std::to_string(e);}

big_int getZero(big_int e){return (big_int)0;}
big_int getOne(big_int e){return (big_int)1;}
