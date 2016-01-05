#include <cgreen/cgreen.h>
#include <Ensemble.h>

static struct BLEnsemble ensemble;

Describe(Ensemble)
BeforeEach(Ensemble) {}
AfterEach(Ensemble) {}

Ensure(Ensemble, initiallyHasZeroSize) {
  blEnsembleInitialize(10, 4, &ensemble);
  assert_that(ensemble.numPtcls, is_equal_to(0));
  blEnsembleFree(&ensemble);
}

Ensure(Ensemble, hasRequestedCapacity) {
  const int requestedCapacity = 6;
  blEnsembleInitialize(requestedCapacity, 4, &ensemble);
  assert_that(ensemble.maxNumPtcls, is_equal_to(requestedCapacity));
  blEnsembleFree(&ensemble);
}

Ensure(Ensemble, canBeFreed) {
  blEnsembleInitialize(5, 4, &ensemble);
  blEnsembleFree(&ensemble);
}

Ensure(Ensemble, createSpace) {
  blEnsembleInitialize(5, 4, &ensemble);
  blEnsembleCreateSpace(2, &ensemble);
  assert_that(ensemble.numPtcls, is_equal_to(2));
  blEnsembleFree(&ensemble);
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Ensemble, initiallyHasZeroSize);
  add_test_with_context(suite, Ensemble, hasRequestedCapacity);
  add_test_with_context(suite, Ensemble, canBeFreed);
  add_test_with_context(suite, Ensemble, createSpace);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}
