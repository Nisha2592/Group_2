cmake_minimum_required(VERSION 3.16)
project(GTESTDEMO VERSION 1.0 LANGUAGES C)

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# Library and executable definition
add_library(mdlib SHARED include/utilities.c include/output.c include/input.c include/force.c include/cleanup.c include/allocation.c)
add_compile_options(-O3 -mfma -ftree-vectorize -fopenmp-simd -ffast-math -fexpensive-optimizations -fomit-frame-pointer -march=native -ffp-contract=fast -msse3 -DLJMD_VERSION=1.0)

add_executable(ljmd src/main.c)
#target_link_libraries(ljmd mdlib m ${MPI_LIBRARIES} OpenMP::OpenMP_C) 
target_link_libraries(ljmd mdlib)

# for MPI
option(ENABLE_MPI "Enable MPI calls" OFF)
if(ENABLE_MPI)
  set(MPI_C_SKIP_MPIC TRUE)
  find_package(MPI REQUIRED)
  if(MPI_FOUND)
    message(STATUS "MPI found.")
    add_compile_definitions(USE_MPI=yes)
    target_link_libraries(mdlib PUBLIC MPI::MPI_C)
  else()
    message(WARNING "MPI not found. Proceeding without MPI support.")
  endif()
endif()

# Add OpenMP support
option(ENABLE_OPENMP "Enable OpenMP calls" OFF)
if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_FOUND)
    message(STATUS "OpenMP found, enabling parallelism.")
    target_link_libraries(mdlib PUBLIC OpenMP::OpenMP_C)
    target_compile_options(mdlib PUBLIC ${OpenMP_C_FLAGS})
  else()
    message(FATAL_ERROR "OpenMP not found. Please install OpenMP support for your compiler.")
  endif()
endif()

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
