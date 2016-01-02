#ifndef __ZELEM_HPP
#define __ZELEM_HPP

#include <string>
#include "types.hpp"

namespace alcp {
    long long add(long long a, long long b, long long p);

    long long russianPeasantMultiplication(long long a, long long b, long long p);

    std::string to_string(big_int e);

    template<typename T>
    bool compatible(T lhs, T rhs);

    template<typename T>
    int unit(T e);

    template<typename T>
    T normalForm(T e);

    template<typename T>
    T getZero(T e);

    template<typename T>
    T getOne(T e);
}
#endif // __ZELEM_HPP
