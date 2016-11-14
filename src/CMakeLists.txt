file(GLOB SRC_FILES *.cpp *.hpp)
add_library(alcp_lib ${SRC_FILES})
add_executable(alcp2 main.cpp)
target_link_libraries (alcp2 alcp_lib)
