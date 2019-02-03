set(HOME_DIR /Users/Owner)

set(ALL_TEST_SRCS
    main_test.cc
    random_number_generator_func_test.cc
    prob_density_func_test.cc
    cumulative_dist_func_test.cc
    black_scholes_explicit_expectation_test.cc
    test_util_test.cc
    euler_maruyama_alg_test.cc
    statistical_functions_test.cc
    )

set(GTEST_INCLUDES
    ${HOME_DIR}/googletest/googletest/include
    ${HOME_DIR}/googletest/googlemock/include
 )
set(LIB_PATHS
    /usr/local/lib
    #path for gtest library
    ${HOME_DIR}/googletest/mybuild/lib
    )


link_directories(${LIB_PATHS})
add_executable(Testset ${ALL_TEST_SRCS})
target_include_directories(Testset PUBLIC ${GTEST_INCLUDES})
target_link_libraries(Testset mylib gtest gtest_main)
