#include <cgreen/cgreen.h>
#include <Utilities.h>
#include <math.h>

Describe(PoissonDist)
BeforeEach(PoissonDist) {}
AfterEach(PoissonDist) {}

Ensure(PoissonDist, isANonNegativeNumber) {
  int n;
  for (n = 0; n < 100; ++n) {
    int i = blGeneratePoisson(0.3);
    assert_that(i, is_greater_than(-1));
  }
}

Ensure(PoissonDist, hasMeanCloseToNBar) {
  int n;
  double nbar = 1.4;
  double mean = 0;
  int numSamples = 1000;
  for (n = 0; n < numSamples; ++n) {
    mean += blGeneratePoisson(nbar);
  }
  mean /= numSamples;
  assert_that_double(fabs(nbar - mean) / nbar, is_less_than_double(0.4));
}

Ensure(PoissonDist, hasVarianceCloseToNBar) {
  int n;
  double nbar = 1.4;
  double variance = 0;
  double mean = 0;
  int numSamples = 10000;
  for (n = 0; n < numSamples; ++n) {
    int sample = blGeneratePoisson(nbar);
    mean += sample;
    variance += sample * sample;
  }
  mean /= numSamples;
  variance /= numSamples;
  variance -= mean * mean;
  assert_that_double(fabs(mean - variance) / mean,
      is_less_than_double(0.4));
}

Ensure(PoissonDist, hasMeanCloseToNBarForLargeN) {
  int n;
  double nbar = 1000.0;
  double mean = 0;
  int numSamples = 1000;
  for (n = 0; n < numSamples; ++n) {
    mean += blGeneratePoisson(nbar);
  }
  mean /= numSamples;
  assert_that_double(fabs(nbar - mean) / nbar,
      is_less_than_double(0.1));
}

Ensure(PoissonDist, hasVarianceCloseToNBarForLargeN) {
  int n;
  double nbar = 1.0e4;
  double variance = 0;
  double mean = 0;
  int numSamples = 1000;
  for (n = 0; n < numSamples; ++n) {
    int sample = blGeneratePoisson(nbar);
    mean += sample;
    variance += sample * sample;
  }
  mean /= numSamples;
  variance /= numSamples;
  variance -= mean * mean;
  printf("mean == %lf\n", mean);
  printf("variance == %lf\n", variance);
  assert_that_double(fabs(mean - variance) / mean,
      is_less_than_double(1.0e-1));
}

int main(int argn, char **argv) {
#ifdef BL_WITH_MPI
  MPI_Init(&argn, &argv);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, PoissonDist, isANonNegativeNumber);
  add_test_with_context(suite, PoissonDist, hasMeanCloseToNBar);
  add_test_with_context(suite, PoissonDist, hasVarianceCloseToNBar);
  add_test_with_context(suite, PoissonDist, hasMeanCloseToNBarForLargeN);
  add_test_with_context(suite, PoissonDist, hasVarianceCloseToNBarForLargeN);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return result;
}

