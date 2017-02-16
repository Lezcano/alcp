#  ALCP [![Build Status](https://travis-ci.org/Lezcano/alcp.svg?branch=master)](https://travis-ci.org/Lezcano/alcp)

> Classical Computer Algebra Algorithms

## Overview

This library gives a modern, efficient and scalable implementation of algebraic data structures and algorithms on them via the use of templates.

The field and ring classes mimic the integral C++ types. This provides an interface that allows the implementation of algorithms that work for both basic integral types and arbitrary rings.

## Building ALCP
In order to build ALCP you will need [CMake][], a C++14 compiler and 
an updated version of the Boost libraries.

In general, you will need to specify a path to a modern C++14 complaint 
compiler using the `-DCMAKE_CXX_COMPILER` variable if the default compiler
is too old.

The library also makes use of the Boost libraries for multiprecision
integers support. In the case that boost libraries are not in the path
you will have to specify the path using the `-DBOOST_ROOT` variable.

Gtest 1.8.0 is already provided in the repository.  Then it is just a 
matter of `cd` to the root directory and setup the build directory
```shell
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/path/to/compiler -DBOOST_ROOT=/path/to/boost
```

The `-DCMAKE_BUILD_TYPE=Release` option will allow to output the highest performance, but it will not check if operations between two elements are compatible (i.e. the elements are in the same field / ring), this might lead into undefined behavior of the library. If we configure the library using just

```shell
cmake ..
```
we will get the version of the library suitable for debugging, that will throw exceptions at runtime if the elements of two operations are not compatible.

We can deactivate the checks passing the argument `-DALCP_NO_CHECKS` to Cmake

## Basic data structures:

Quotient of an ED by a principal ideal: R / \<a\>

- Finite field GF(p) = Z / \<p\>

- Finite field GF(q) = GF(p)[X] / \<f\>

Polynomial ring of an ED: R[X]

- Polynomial ring GF(p)[X]

- Polynomial ring GF(q)[X]

- Polynomial ring Z[X]

## Main algorithms:

- Extended euclides algorithm for an ED

- Chinese remainder algorithm

- Berlekamp

- Cantor-Zassenhaus

- Hensel

- Modular Euclidean algorithm

- BCH codification. General Berlekamp-Massey algorithm. BCH decodification using Berlekamp-Massey
## Misc:

- Addition and multiplication of 63-bit numbers using 64-bit registers - Russian peasant algorithm

- Irreducibility criterion for GF(p)[X]

- Pollard's ρ algorithm for integer factorization

- Discrete logarithm in GF(q) - Pollard's ρ algorithm for logarithms

- Miller-Rabin

# Developers:
This project was developed by [Lezcano](https://github.com/Lezcano) and [damaru2](https://github.com/damaru2)
<!-- Links -->
[CMake]: http://www.cmake.org
