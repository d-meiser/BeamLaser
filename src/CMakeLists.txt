include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${SPRNG_INCLUDE_DIR}
    )

set(BeamLaserLibSources
    Ensemble.c
    Integrator.c
    Partition.c
   )

add_library(beamlaserlib ${BeamLaserLibSources})

