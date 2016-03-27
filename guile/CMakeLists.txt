find_package(Guile)

include_directories(
    ${PROJECT_SOURCE_DIR}/include 
    ${PROJECT_BINARY_DIR}
    ${GUILE_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_PATH}
    )

add_library(guile_beam_laser SHARED guile_beam_laser.c)
target_link_libraries(guile_beam_laser beamlaserlib ${GUILE_LIBRARIES})
install(TARGETS guile_beam_laser DESTINATION lib)
