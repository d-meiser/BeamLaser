include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    )
if (BL_WITH_MPI)
  include_directories(${MPI_C_INCLUDE_PATH} SYSTEM)
endif ()
add_executable(BeamLaser BeamLaser.c)
target_link_libraries(BeamLaser beamlaserlib ${SPRNG_LIB} m)
if (BL_WITH_MPI)
  target_link_libraries(BeamLaser ${MPI_C_LIBRARIES})
endif ()

add_test(BeamLaser BeamLaser)

