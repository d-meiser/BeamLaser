find_package(Guile)

include_directories(
    ${PROJECT_SOURCE_DIR}/include 
    ${PROJECT_BINARY_DIR}
    )
include_directories(${GUILE_INCLUDE_DIRS} SYSTEM)
if(BL_WITH_MPI)
  include_directories(${MPI_C_INCLUDE_PATH} SYSTEM)
endif()

add_library(guile_beam_laser SHARED guile_beam_laser.c)
target_link_libraries(guile_beam_laser beamlaserlib ${GUILE_LIBRARIES})
install(TARGETS guile_beam_laser DESTINATION lib)
