include_directories(
    ${PROJECT_SOURCE_DIR}/include
    )
add_executable(BeamLaser BeamLaser.c)
target_link_libraries(BeamLaser beamlaserlib ${SPRNG_LIB} m)

add_test(BeamLaser BeamLaser)

