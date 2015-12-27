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

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Integrator,
    yieldsErrorWhenCreatedFromUnknownName);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

