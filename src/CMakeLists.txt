include_directories(${PROJECT_SOURCE_DIR}/include)

set(BeamLaserLibSources
    Ensemble.c
    RingBuffer.c
    Partition.c
   )

add_library(beamlaserlib ${BeamLaserLibSources})

add_executable(beamlaser beamlaser.c)
