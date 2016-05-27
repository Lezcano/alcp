# Mostyl based on the CMake from Manu343726
# https://github.com/GueimUCM/siplasplas/blob/master/CMakeLists.txt

cmake_minimum_required(VERSION 2.8.8)
set(PROJECT_NAME_STR alcp)
project(${PROJECT_NAME_STR})

find_package(Threads REQUIRED)

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
    set(CPP_GLOBAL_DEBUG_COMPILE_OPTIONS -O0 -g3)
    set(CPP_GLOABL_RELEASE_COMPILE_OPTIONS -O3 -g0 -DALCP_NO_CHECKS)
else()
    set(CPP_GLOBAL_ALL_COMPILE_OPTIONS /DCTTI_STRING_MAX_LENGTH=512)
endif()

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)

include_directories(${COMMON_INCLUDES})
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

#-------------------
# Test
#-------------------
enable_testing()
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)
file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/tests/*.cpp)
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_test(test1 ${PROJECT_TEST_NAME})

target_link_libraries(${PROJECT_TEST_NAME} ${CMAKE_THREAD_LIBS_INIT})

add_executable(${PROJECT_NAME_STR} ${SRC_FILES})
