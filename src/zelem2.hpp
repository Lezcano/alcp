#ifndef __ZELEM_HPP
#define __ZELEM_HPP

#include <string>
#include <sstream>

#include "types.hpp"

namespace alcp {
    long long add(long long a, long long b, long long p);

    long long russianPeasantMultiplication(long long a, long long b, long long p);

    // Enable if it is multiprecision
    template<class Int>
    typename std::enable_if<boost::multiprecision::is_number<Int>::value, std::string>::type
    to_string (Int e)
    {
        return e.template convert_to<std::string>();
    }

    // Enable if it is basic type
    template<class Int>
    typename std::enable_if<std::is_integral<Int>::value, std::string>::type
    to_string(Int e)
    {
        return std::to_string(e);
    }

    template<typename T, class = typename std::enable_if<is_integral<T>::value>::type>
    std::string to_string_coef(T e){
        if(e < 0) return to_string(e);
        return "+" + to_string(e);
    }

    template<typename T, class = typename std::enable_if<is_integral<T>::value>::type>
    bool compatible(T, T) {
        return true;
    }

    template<typename T, class = typename std::enable_if<is_integral<T>::value>::type>
    T unit(T e) {
        return e >= 0 ? 1 : -1;
    }

    template<typename T, class = typename std::enable_if<is_integral<T>::value>::type>
    T normalForm(T e) {
        return e / unit<T>(e);
    }

    template<typename T, class = typename std::enable_if<is_integral<T>::value>::type>
    T getZero(T) {
        return 0;
    }

    template<typename T, class = typename std::enable_if<is_integral<T>::value>::type>
    T getOne(T) {
        return 1;
    }

}
#endif // __ZELEM_HPP
