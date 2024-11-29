cmake_minimum_required(VERSION 3.16)
project(GTESTDEMO VERSION 1.0 LANGUAGES C)

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

add_library(mdlib ../include/utilities.c ../include/force.c)

option(ENABLE_TESTING "Enable building unit tests" OFF)
if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)

  add_executable(test_force test_force.cpp)
  target_include_directories(test_force PRIVATE ${CMAKE_SOURCE_DIR}) 
  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME force COMMAND test_force)

  add_executable(test_verlet test_verlet.cpp)
  target_include_directories(test_verlet PRIVATE ${CMAKE_SOURCE_DIR})
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Verlet COMMAND test_verlet)

  add_executable(test_utils test_utilities.cpp)
  target_include_directories(test_utils PRIVATE ${CMAKE_SOURCE_DIR})
  target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Utils COMMAND test_utils)
endif()
