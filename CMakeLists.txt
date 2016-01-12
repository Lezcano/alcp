cmake_minimum_required(VERSION 3.3)
project(alcp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES
    src/exceptions.hpp
    src/factorizationFq.hpp
    src/fpelem.hpp
    src/fpxelem.hpp
    src/fqelem.hpp
    src/fqxelem.hpp
    src/generalPurpose.cpp
    src/generalPurpose.hpp
    src/hensel.cpp
    src/hensel.hpp
    src/henselSubsets.cpp
    src/henselSubsets.hpp
    src/integerCRA.cpp
    src/integerCRA.hpp
    src/main.cpp
    src/modularGCD.cpp
    src/modularGCD.hpp
    src/polRing.hpp
    src/quotientRing.hpp
    src/types.hpp
    src/userInterface.hpp
    src/userInterface.cpp
    src/zelem.cpp
    src/zelem.hpp
    src/zxelem.hpp)

add_executable(alcp ${SOURCE_FILES})
