set(ALL_SRCS
    true_value.cc)

#Make sure that the dependencies among headers are in a correct order.
set(ALL_INCLUDES
    ${PROJECT_SOURCE_DIR}/include
    /usr/local/include
    /usr/local/include/boost
    )

set(LIB_PATHS
    /usr/local/lib)

add_library(mylib STATIC ${ALL_SRCS})
target_include_directories(mylib PUBLIC ${ALL_INCLUDES})
link_directories(${LIB_PATHS})
target_link_libraries(mylib QuantLib)
