include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}
    ${CGREEN_PUBLIC_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_PATH}
    )

set(TESTS
    testBSP
    testDipoleOperatorTLA
    testEnsemble
    testIntegrator
    testParticleSource
    testPartition
    testPushParticles
    testRemoveParticles
    )

foreach(t ${TESTS})
  add_executable(${t} ${t}.c)
  target_link_libraries(${t} beamlaserlib cgreen_shared ${SPRNG_LIB})
  if(BL_WITH_MPI)
    target_link_libraries(${t} ${MPI_C_LIBRARIES})
  endif()
  add_test(${t} ${t})
endforeach()


if(BL_WITH_MPI)
  set(MPI_TESTS
      testParticleSourceMPI
      testBcast
      testAddAll
      )

  foreach(t ${MPI_TESTS})
    add_executable(${t} ${t}.c)
    target_link_libraries(${t}
                          beamlaserlib
                          cgreen_shared
                          ${SPRNG_LIB}
                          ${MPI_C_LIBRARIES}
                         )
    add_test(${t}_1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ./${t})
    add_test(${t}_2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ./${t})
  endforeach()

endif()
