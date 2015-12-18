cmake_minimum_required(VERSION 2.6.2)
if (NOT(CMAKE_MAJOR_VERSION LESS 3))
  cmake_policy(SET CMP0048 OLD)
endif ()
project(BeamLaser LANGUAGES C)

enable_testing()

option(BL_BUILD_TESTS "Whether to build tests" ON)

add_subdirectory(src)
if(BL_BUILD_TESTS)
  add_subdirectory(cgreen)
  add_subdirectory(tests)
endif()

