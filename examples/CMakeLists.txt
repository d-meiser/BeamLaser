add_executable(BeamLaser BeamLaser.c)
target_link_libraries(BeamLaser beamlaserlib ${SPRNG_LIB})

add_test(BeamLaser BeamLaser)

