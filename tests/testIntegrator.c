#include <cgreen/cgreen.h>
#include <Integrator.h>
#include <math.h>

Describe(Integrator);
BeforeEach(Integrator) {
  significant_figures_for_assert_double_are(6);
}
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

struct DecayCtx {
  double gamma;
};

void exponentialDecay(double t, int n, const double* x,
                      double* y, void* ctx) {
  (void)t;
  struct DecayCtx* decCtx = (struct DecayCtx*)ctx;
  int i;
  for (i = 0; i < n; ++i) {
    y[i] = -decCtx->gamma * x[i];
  }
}

Ensure(Integrator, accuratelyIntegratesExponentialDecay) {
  BLIntegrator integrator;
  blIntegratorCreate("RK4", 1, &integrator);
  double x = 3.7;
  double y = 0.0;
  struct DecayCtx ctx = {2.0};
  double dt = 0.001;
  blIntegratorTakeStep(integrator, 0.0, dt, exponentialDecay, &x, &y, &ctx);
  assert_that_double(y, is_equal_to_double(x * exp(-dt * ctx.gamma)));
  blIntegratorDestroy(&integrator);
}

Ensure(Integrator, worksForVectors) {
  BLIntegrator integrator;
  blIntegratorCreate("RK4", 2, &integrator);
  double x[2] = {3.7, 1.2};
  double y[2] = {0.0, 0.0};
  struct DecayCtx ctx = {2.0};
  double dt = 0.001;
  blIntegratorTakeStep(integrator, 0.0, dt, exponentialDecay, x, y, &ctx);
  assert_that_double(y[0], is_equal_to_double(x[0] * exp(-dt * ctx.gamma)));
  assert_that_double(y[1], is_equal_to_double(x[1] * exp(-dt * ctx.gamma)));
  blIntegratorDestroy(&integrator);
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Integrator,
    yieldsErrorWhenCreatedFromUnknownName);
  add_test_with_context(suite, Integrator, canBeCreateFromKnownName);
  add_test_with_context(suite, Integrator, accuratelyIntegratesExponentialDecay);
  add_test_with_context(suite, Integrator, worksForVectors);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

