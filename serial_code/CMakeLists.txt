cmake_minimum_required(VERSION 3.16)
project(GTESTDEMO VERSION 1.0 LANGUAGES C)

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# Library and executable definition
add_library(mdlib SHARED include/utilities.c include/output.c include/input.c include/force.c include/cleanup.c include/allocation.c)
target_link_libraries(mdlib PUBLIC m) 
add_compile_options(-O3 -mfma -ftree-vectorize -fopenmp-simd -ffast-math -fexpensive-optimizations -fomit-frame-pointer -march=native -ffp-contract=fast -msse3 -DLJMD_VERSION=1.0)

add_executable(ljmd src/main.c)
target_link_libraries(ljmd m mdlib)
