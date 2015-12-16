include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(beamlaserlib ensemble.c)

add_executable(beamlaser beamlaser.c)
