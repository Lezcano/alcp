#ifndef __TYPES_HPP
#define __TYPES_HPP

#include <boost/multiprecision/cpp_int.hpp>
#include <type_traits>

namespace alcp {
    //using big_int = long long int;
    //using big_int = boost::multiprecision::cpp_int;
    using big_int =  boost::multiprecision::int256_t;



    template<class T> struct is_integral :
        std::integral_constant < bool,
            std::is_integral<T>::value ||
            std::is_same<typename std::decay<T>::type, big_int>::value
        >{};
}

#endif // __TYPES_HPP
