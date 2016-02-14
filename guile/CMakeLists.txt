find_package(Guile)

include_directories(
    ${PROJECT_SOURCE_DIR}/include 
    ${PROJECT_BINARY_DIR}
    ${GUILE_INCLUDE_DIRS})

add_library(bl_guile SHARED bl_guile.c)
target_link_libraries(bl_guile beamlaserlib ${GUILE_LIBRARIES})
