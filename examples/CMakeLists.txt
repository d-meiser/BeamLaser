include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    ${MPI_C_INCLUDE_PATH}
    )
if (BL_WITH_MPI)
  include_directories(${MPI_C_INCLUDE_PATH} SYSTEM)
endif ()
add_executable(BeamLaserTLA BeamLaserTLA.c)
target_link_libraries(BeamLaserTLA beamlaserlib ${SPRNG_LIB} m)
if (BL_WITH_MPI)
  target_link_libraries(BeamLaserTLA ${MPI_C_LIBRARIES})
endif ()

add_test(BeamLaserTLA BeamLaserTLA)

