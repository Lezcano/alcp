#ifndef __ZELEM_HPP
#define __ZELEM_HPP

#include <string>
#include <sstream>

#include "types.hpp"

namespace alcp {
    long long add(long long a, long long b, long long p);

    long long russianPeasantMultiplication(long long a, long long b, long long p);

    // Enable if it is boost multiprecision
    template<class Int>
    typename std::enable_if<bmp::is_number<Int>::value, std::string>::type
    to_string (const Int& e)
    {
        return e.template convert_to<std::string>();
    }

    // Enable if is base type
    template<class Int>
    std::enable_if_t<std::is_integral<Int>::value, std::string>
    to_string(const Int& e)
    {
        return std::to_string(e);
    }

    template<class Int>
    std::string to_string_coef(const Int& e)
    {
        static_assert(is_integral<Int>::value, "Type is not a supported integer.");
        return e < 0 ? to_string(e) : "+" + to_string(e);
    }

    template<class Int>
    Int compatible(const Int&, const Int&) {
        static_assert(is_integral<Int>::value, "Type is not a supported integer.");
        return true;
    }

    template<class Int>
    Int unit(const Int& e) {
        static_assert(is_integral<Int>::value, "Type is not a supported integer.");
        return e >= 0 ? 1 : -1;
    }

    template<class Int>
    Int normalForm(const Int& e) {
        static_assert(is_integral<Int>::value, "Type is not a supported integer.");
        return e / unit<Int>(e);
    }

    template<class Int>
    Int getZero(Int) {
        static_assert(is_integral<Int>::value, "Type is not a supported integer.");
        return 0;
    }

    template<class Int>
    Int getOne(Int) {
        static_assert(is_integral<Int>::value, "Type is not a supported integer.");
        return 1;
    }

}
#endif // __ZELEM_HPP
