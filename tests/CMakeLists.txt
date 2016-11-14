cmake_minimum_required(VERSION 3.3)

set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)

find_package(GTest REQUIRED)
include_directories(${GTEST_SOURCE_DIR}/include ${GTEST_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src)

add_executable(${PROJECT_TEST_NAME} tests.cpp)

target_link_libraries(${PROJECT_TEST_NAME}
    gtest
    gtest_main
    ${CMAKE_THREAD_LIBS_INIT}
    alcp_lib)

enable_testing()
add_test(test1 ${PROJECT_TEST_NAME})
