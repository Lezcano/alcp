#include <string>

#include "zelem.hpp"
// Types inlcuded in zelem.hpp

/**
 * Addition of 63 bit numbers without overflow
 */
long long add(long long a, long long b, long long p){
    return static_cast<long long>(
            static_cast<unsigned long long>(a) +
            static_cast<unsigned long long>(b) %
            static_cast<unsigned long long>(p));
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
        if (a & 1)
            res = (res + static_cast<unsigned long long>(b)) % static_cast<unsigned long long>(p);
        a >>= 1;
        b = (b << 1) % p;
    }
    return static_cast<long long>(res);
}

//std::string to_string(big_int e){return e.convert_to<std::string>();}

template <typename T>
bool compatible(T lhs,T rhs){
    (void) lhs; (void) rhs; // Silences -Wunused-param and does nothng
    return true;
}

template <typename T>
int unit(T e){
    return e >= 0 ? 1 : -1;
}

template <typename T>
T normalForm(T e){ return e/unit<T>(e);}

template <typename T>
T getZero(T e){
    (void) e; // Silences -Wunused-param and does nothng
    return 0;
}
template <typename T>
T getOne(T e){
    (void) e; // Silences -Wunused-param and does nothng
    return 1;
}


template bool compatible<big_int>(big_int lhs, big_int rhs);
template int unit<big_int>(big_int e);
template big_int normalForm<big_int>(big_int e);
template big_int getZero<big_int>(big_int e);
template big_int getOne<big_int>(big_int e);

template bool compatible<int>(int lhs, int rhs);
template int unit<int>(int e);
template int normalForm<int>(int e);
template int getZero<int>(int e);
template int getOne<int>(int e);

//template bool compatible<long long>(long long lhs, long long rhs);
//template long long unit<long long>(long long e);
//template long long normalForm<long long>(long long e);
//template long long getZero<long long>(long long e);
//template long long getOne<long long>(long long e);


