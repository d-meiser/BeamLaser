include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CGREEN_PUBLIC_INCLUDE_DIRS}
    )

set(TESTS
    testRingBuffer)

foreach(t ${TESTS})
  add_executable(${t} ${t}.c)
  target_link_libraries(${t} beamlaserlib cgreen_shared)
  add_test(${t} ${t})
endforeach()

