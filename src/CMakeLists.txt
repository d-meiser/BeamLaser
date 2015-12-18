include_directories(${PROJECT_SOURCE_DIR}/include)

set(BeamLaserLibSources
    ensemble.c
    RingBuffer.c
    Partition.c
   )

add_library(beamlaserlib ${BeamLaserLibSources})

add_executable(beamlaser beamlaser.c)
