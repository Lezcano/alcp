set(PROJECT_TEST_NAME ${PROJECT_NAME}_test)

find_package(Threads REQUIRED)

add_executable(${PROJECT_TEST_NAME} tests.cpp)

include_directories(${PROJECT_SOURCE_DIR}/src)

target_link_libraries(${PROJECT_TEST_NAME}
    gtest
    gtest_main
    ${CMAKE_THREAD_LIBS_INIT}
    ${PROJECT_LIB})

add_test(all_tests ${PROJECT_TEST_NAME})
