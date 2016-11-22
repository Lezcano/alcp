# Mostyl based on the CMake from Manu343726
# https://github.com/GueimUCM/siplasplas/blob/master/CMakeLists.txt

cmake_minimum_required(VERSION 2.8.6)

set(PROJECT_NAME alcp)
set(PROJECT_LIB ${PROJECT_NAME}_lib)


project(${PROJECT_NAME})


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Project options
if(NOT MSVC)
    # gtest needs fileno() function at file gtest/include/gtest/internal/gtest-port.h,
    # but is considered an extension in cygwin GCC setup. Use GNU extensions.
    # See http://stackoverflow.com/questions/18784112/googletest-1-6-with-cygwin-1-7-compile-error-fileno-was-not-declared-in-this
    if(CYGWIN)
        set(STD_CXX gnu++14)
    else()
        set(STD_CXX c++14)
    endif()

    set(CPP_GLOBAL_ALL_COMPILE_OPTIONS -std=${STD_CXX} -Wall -pedantic )
    set(CPP_GLOBAL_DEBUG_COMPILE_OPTIONS -O0 -g3 -Wextra -Werror -fsanitize=undefined -fsanitize=address -fsanitize=thread -fsanitize=memory)
    set(CPP_GLOABL_RELEASE_COMPILE_OPTIONS -O3 -g0 -DALCP_NO_CHECKS)
else()
    set(CPP_GLOBAL_ALL_COMPILE_OPTIONS /DCTTI_STRING_MAX_LENGTH=512)
endif()

add_subdirectory(src)

# Testing
add_subdirectory(deps/gtest-1.8.0)
include(CTest)
add_subdirectory(tests)
