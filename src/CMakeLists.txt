cmake_minimum_required(VERSION 3.3)
project(src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    exceptions.hpp
    factorizationFq.cpp
    factorizationFq.hpp
    fp.cpp
    fp.hpp
    fpelem.cpp
    fpelem.hpp
    fpxelem.cpp
    fpxelem.hpp
    fq.cpp
    fq.hpp
    fqelem.cpp
    fqelem.hpp
    fqxelem.cpp
    fqxelem.hpp
    generalPurpose.cpp
    generalPurpose.hpp
    hensel.cpp
    hensel.hpp
    henselSubsets.cpp
    henselSubsets.hpp
    integerCRA.cpp
    integerCRA.hpp
    main.cpp
    #modularGCD.cpp
    polRing.hpp
    types.hpp
    zelem.cpp
    zelem.hpp
    zxelem.cpp
    zxelem.hpp)

add_executable(src ${SOURCE_FILES})