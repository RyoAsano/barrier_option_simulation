set(HOME_DIR /Users/Owner)

set(ALL_MAIN_SRCS
    black_scholes_upper_barrier_main.cc
   )

set(LIB_PATHS
    /usr/local/lib
    )


link_directories(${LIB_PATHS})
add_executable(black_scholes_upper_barrier_simulation ${ALL_MAIN_SRCS})
target_link_libraries(black_scholes_upper_barrier_simulation mylib)
