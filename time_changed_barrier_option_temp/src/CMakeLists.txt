set(ALL_SRCS 
    function.cpp
    generator_Euler_Maruyama.cpp
    generator_Euler_Maruyama_conditional.cpp
    generator_Euler_Maruyama_merger.cpp
    generator_first_hitting_time.cpp
    generator_rand_num.cpp
    vector_field_class.cpp
    )

set(ALL_INCLUDES
    ${PROJECT_SOURCE_DIR}/include
    /usr/local/include
    /usr/local/include/boost
    )

set(ALL_LIB
    /usr/local/lib
    )

link_directories(${ALL_LIB})
add_library(mylib STATIC ${ALL_SRCS})
target_include_directories(mylib PUBLIC ${ALL_INCLUDES})
target_link_libraries(mylib QuantLib)
