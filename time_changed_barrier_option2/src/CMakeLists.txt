set(ALL_SRCS
    cumulative_dist_func.cc
    euler_maruyama_alg.cc
    prob_density_func.cc
    random_number_generator_func.cc
    black_scholes_explicit_expectation.cc
    test_util.cc
    )

#Make sure that the dependencies among headers are in a correct order.
set(ALL_INCLUDES
    ${PROJECT_SOURCE_DIR}/include
    /usr/local/include
    /usr/local/include/boost
   )

set(LIB_PATHS
    /usr/local/lib
    )

link_directories(${LIB_PATHS})
add_library(mylib STATIC ${ALL_SRCS})
target_include_directories(mylib PUBLIC ${ALL_INCLUDES})
target_link_libraries(mylib QuantLib)
