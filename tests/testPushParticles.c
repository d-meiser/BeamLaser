#include <cgreen/cgreen.h>
#include <Ensemble.h>

static struct BLEnsemble ensemble;

Describe(Push)
BeforeEach(Push) {
  static const int numPtcls = 5;
  blEnsembleInitialize(numPtcls + 1, 2, &ensemble);
  int i;
  for (i = 0; i < numPtcls; ++i) {
    ensemble.x[i] = 0;
    ensemble.y[i] = 0;
    ensemble.z[i] = 0;
    ensemble.vx[i] = rand() / (double)RAND_MAX;
    ensemble.vy[i] = rand() / (double)RAND_MAX;
    ensemble.vz[i] = rand() / (double)RAND_MAX;
  }
  ensemble.numPtcls = numPtcls;
}
AfterEach(Push) {
  blEnsembleFree(&ensemble);
}

Ensure(Push, movesParticlesTheRightDistance) {
  blEnsemblePush(0.7, &ensemble);
  int i;
  for (i = 0; i != ensemble.numPtcls; ++i) {
    assert_that_double(ensemble.x[i], is_equal_to_double(ensemble.vx[i] * 0.7));
    assert_that_double(ensemble.y[i], is_equal_to_double(ensemble.vy[i] * 0.7));
    assert_that_double(ensemble.z[i], is_equal_to_double(ensemble.vz[i] * 0.7));
  }
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Push, movesParticlesTheRightDistance);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

