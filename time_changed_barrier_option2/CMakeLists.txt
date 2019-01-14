cmake_minimum_required(VERSION 2.8)
project(time_changed_barrier_option CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_definitions(-std=c++11)

add_subdirectory(src)
add_subdirectory(test)
