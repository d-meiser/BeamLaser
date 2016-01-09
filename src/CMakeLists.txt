include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    ${SPRNG_INCLUDE_DIR}
    )

set(BeamLaserLibSources
    DipoleOperator.c
    Ensemble.c
    Integrator.c
    ParticleSource.c
    Partition.c
    Utilities.c
   )

add_library(beamlaserlib ${BeamLaserLibSources})

