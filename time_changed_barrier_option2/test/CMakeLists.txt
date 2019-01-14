set(ALL_LIB
    /usr/local/lib
    )

link_directories(${ALL_LIB})

add_executable(Main main.cc)
target_link_libraries(Main mylib)
