cmake_minimum_required(VERSION 3.13.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_definitions(-std=c++11)

project(test CXX)
add_subdirectory(src)
add_subdirectory(test)
