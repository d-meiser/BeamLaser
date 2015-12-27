#include <cgreen/cgreen.h>
#include <Integrator.h>

Describe(Integrator);
BeforeEach(Integrator) {}
AfterEach(Integrator) {}

Ensure(Integrator, yieldsErrorWhenCreatedFromUnknownName) {
  BLIntegrator integrator;
  BL_STATUS stat = blIntegratorCreate("Unknown name", 2, &integrator);
  assert_that(stat, is_equal_to(BL_UNKNOWN_INTEGRATOR));
}

Ensure(Integrator, canBeCreateFromKnownName) {
  BLIntegrator integrator;
  BL_STATUS stat = blIntegratorCreate("RK4", 2, &integrator);
  assert_that(stat, is_equal_to(BL_SUCCESS));
  assert_that(integrator, is_not_null);
  blIntegratorDestroy(&integrator);
}
int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Integrator,
    yieldsErrorWhenCreatedFromUnknownName);
  add_test_with_context(suite, Integrator, canBeCreateFromKnownName);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

