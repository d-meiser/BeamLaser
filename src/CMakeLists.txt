include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    ${SPRNG_INCLUDE_DIR}
    )
if(BL_WITH_MPI)
  include_directories(${MPI_C_INCLUDE_PATH} SYSTEM)
endif()

set(BeamLaserLibSources
    AtomFieldInteraction.c
    Diagnostics.c
    DipoleOperator.c
    Ensemble.c
    Integrator.c
    ModeFunction.c
    ParticleSource.c
    Partition.c
    Update.c
    Utilities.c
   )

add_library(beamlaserlib ${BeamLaserLibSources})

