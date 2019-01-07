cmake_minimum_required(VERSION 2.8)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_definitions(-std=c++11)

project(test CXX)
include_directories(/usr/local/include)
include_directories(/usr/local/include/boost)
link_directories(/usr/local/lib)
add_executable(Main
	main.cpp
	BlackScholesTruePrice.cpp
	generator_first_hitting_time.cpp
	function.cpp
	generator_rand_num.cpp
	generator_Euler_Maruyama.cpp
	main.cpp
	generator_Euler_Maruyama_conditional.cpp
	generator_Euler_Maruyama_merger.cpp
	vector_field_class.cpp)
target_link_libraries(Main QUantLib)
