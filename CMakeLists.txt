cmake_minimum_required(VERSION 3.16)
project(GTESTDEMO VERSION 1.0 LANGUAGES C)

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# for MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_library(mdlib SHARED include/utilities.c include/output.c include/input.c include/force.c include/cleanup.c include/allocation.c)

add_compile_options(-Wall -O3 -ffast-math -fexpensive-optimizations -ffp-contract=fast -msse3 -DLJMD_VERSION=1.0)
# add_compile_options(-DLJMD_VERSION=1.0 -pg)
add_executable(ljmd src/main.c)

target_link_libraries(ljmd mdlib m ${MPI_LIBRARIES}) 

option(ENABLE_TESTING "Enable building unit tests" OFF)

if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)
  add_executable(test_verlet tests/test_verlet.cpp)
  target_include_directories(test_verlet PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Verlet COMMAND test_verlet)

  add_executable(test_utils tests/test_utils.cpp)
  target_include_directories(test_utils PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Utils COMMAND test_utils)
endif()
