include_directories(
    ${PROJECT_SOURCE_DIR}/include 
    ${PROJECT_BINARY_DIR}
    )
if(BL_WITH_MPI)
  include_directories(${MPI_C_INCLUDE_PATH} SYSTEM)
endif()
find_package(Guile)
include_directories(${GUILE_INCLUDE_DIRS} SYSTEM)

add_library(guile_beam_laser SHARED guile_beam_laser.c)
target_link_libraries(guile_beam_laser beamlaserlib ${GUILE_LIBRARIES})
install(TARGETS guile_beam_laser DESTINATION lib)

if(BL_BUILD_TESTS)
  add_subdirectory(tests)
endif()
