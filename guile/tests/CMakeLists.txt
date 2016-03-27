add_test(
    NAME test_guile_ensemble
    COMMAND guile -L ${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_CURRENT_SOURCE_DIR}/test-ensemble.scm)

