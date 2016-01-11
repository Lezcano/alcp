#include "zelem.hpp"


namespace alcp {
/**
 * Addition of 63 bit numbers without overflow
 */
    long long add(long long a, long long b, long long p) {
        return (static_cast<unsigned long long>(a)
                + static_cast<unsigned long long>(b))
               % static_cast<unsigned long long>(p);
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
    long long russianPeasantMultiplication(long long a, long long b, long long p) {
        unsigned long long res = 0;
        while (a != 0) {
            if (a & 1)
                res = (res + static_cast<unsigned long long>(b))
                      % static_cast<unsigned long long>(p);
            a >>= 1;
            b = (b << 1) % p;
        }
        return static_cast<long long>(res);
    }

    std::string to_string(big_int e){return e.convert_to<std::string>();}
}


