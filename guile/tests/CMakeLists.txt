file(COPY test-ensemble.scm DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(
    NAME test_guile_ensemble
    COMMAND guile test-ensemble.scm)

