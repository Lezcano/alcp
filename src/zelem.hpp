#ifndef __ZELEM_HPP
#define __ZELEM_HPP

#include <string>
#include <sstream>

#include "types.hpp"

namespace alcp {
    long long add(long long a, long long b, long long p);

    long long russianPeasantMultiplication(long long a, long long b, long long p);

    // Enable if is multiprecision
    template<class Int>
    typename std::enable_if<boost::multiprecision::is_number<Int>::value, std::string>::type
    to_string (const Int& e)
    {
        return e.template convert_to<std::string>();
    }

    // Enable if is base type
    template<class Int>
    typename std::enable_if<std::is_integral<Int>::value, std::string>::type
    to_string(const Int& e)
    {
        return std::to_string(e);
    }

    template<class Int>
    typename std::enable_if<is_integral<Int>::value, std::string>::type
    to_string_coef(const Int& e){
        if(e < 0) return to_string(e);
        return "+" + to_string(e);
    }

    template<class Int>
    typename std::enable_if<is_integral<Int>::value, Int>::type
    compatible(const Int&, const Int&) {
        return true;
    }

    template<class Int>
    typename std::enable_if<is_integral<Int>::value, Int>::type
    unit(const Int& e) {
        return e >= 0 ? 1 : -1;
    }

    template<class Int>
    typename std::enable_if<is_integral<Int>::value, Int>::type
    normalForm(const Int& e) {
        return e / unit<Int>(e);
    }

    template<class Int>
    typename std::enable_if<is_integral<Int>::value, Int>::type
    getZero(Int) {
        return 0;
    }

    template<class Int>
    typename std::enable_if<is_integral<Int>::value, Int>::type
    getOne(Int) {
        return 1;
    }

}
#endif // __ZELEM_HPP
