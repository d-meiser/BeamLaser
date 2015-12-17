cmake_minimum_required(VERSION 2.6.2)
cmake_policy(SET CMP0048 NEW)
project(BeamLaser VERSION 0.0.0.0)

option(BL_BUILD_TESTS "Whether to build tests" ON)

if(BL_BUILD_TESTS)
  add_subdirectory(cgreen)
endif()

add_subdirectory(src)
