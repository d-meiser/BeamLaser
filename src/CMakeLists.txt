include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${SPRNG_INCLUDE_DIR}
    )

set(BeamLaserLibSources
    Ensemble.c
    Integrator.c
    ParticleSource.c
    Partition.c
    Utilities.c
   )

add_library(beamlaserlib ${BeamLaserLibSources})

