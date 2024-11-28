cmake_minimum_required(VERSION 3.16)
project(GTESTDEMO VERSION 1.0 LANGUAGES C)

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# Library and executable definition
add_library(mdlib SHARED include/utilities.c include/output.c include/input.c include/force.c include/cleanup.c include/allocation.c)
add_compile_options(-Wall -DOPENMP -O3 -mfma -ftree-vectorize -fopenmp-simd -ffast-math -fexpensive-optimizations -fomit-frame-pointer -march=native  -ffp-contract=fast -msse3 -DLJMD_VERSION=1.0)

add_executable(ljmd src/main.c)
#target_link_libraries(ljmd mdlib m ${MPI_LIBRARIES} OpenMP::OpenMP_C) 
target_link_libraries(ljmd mdlib)

# for MPI
option(ENABlE_MPI "Enable MPI calls" OFF)
if(ENABLE_MPI)
 set(MPI_C_SKIP_MPIC TRUE)
 find_package(MPI REQUIRED)
 if (MPI_FOUND)
  message(STATUS "MPI found.")
  add_definitions(-DMPI_ENABLED)
  include_directories(${MPI_INCLUDE_PATH})
  set(LIBS ${LIBS} MPI::MPI_C)
  target_link_libraries(mdlib PUBLIC MPI::MPI_C)
  #set(LIBS ${LIBS} MPI::MPI_C)  # Add MPI library
  endif()
else()
  message(WARNING "MPI not found. Proceeding without MPI support.")
endif()

# Add OpenMP support
option(ENABLE_OPENMP "Enable OpenMP calls" OFF)
if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_FOUND)
    message(STATUS "OpenMP found, enabling parallelism")
    target_link_libraries(mdlib PUBLIC OpenMP::OpenMP_C)
    #add_compile_options(${OpenMP_C_FLAGS})
    #add_link_options(${OpenMP_C_FLAGS})
  endif()
else()
  message(FATAL_ERROR "OpenMP not found. Please install OpenMP support for your compiler.")
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
