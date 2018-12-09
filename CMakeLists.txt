cmake_minimum_required(VERSION 2.6)
project(jacobi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
add_executable(jacobi main.cpp)

install(TARGETS jacobi RUNTIME DESTINATION bin)
