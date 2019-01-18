set(HOME_DIR /Users/Owner)

set(ALL_TEST_SRCS
    main_test.cc
    true_value_test.cc
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
