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
    FieldUpdate.c
    Integrator.c
    ModeFunction.c
    ParticleSource.c
    Partition.c
    PushUpdate.c
    Sink.c
    Source.c
    Update.c
    Utilities.c
   )

message(STATUS "CMAKE_POSITION_INDEPENDENT_CODE == ${CMAKE_POSITION_INDEPENDENT_CODE}")
add_library(beamlaserlib ${BeamLaserLibSources})

