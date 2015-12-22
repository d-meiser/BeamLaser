#include <cgreen/cgreen.h>
#include <Ensemble.h>

static struct BLEnsemble ensemble;

Describe(RemoveBelow)
BeforeEach(RemoveBelow) {
  static const int numPtcls = 20;
  blEnsembleInitialize(numPtcls + 1, 2, &ensemble);
  int i;
  for (i = 0; i < numPtcls; ++i) {
    ensemble.x[i] = rand() / (double)RAND_MAX;
  }
  ensemble.buffer.end = numPtcls;
}
AfterEach(RemoveBelow) {
  blEnsembleFree(&ensemble);
}

Ensure(RemoveBelow, doesNotIncreaseNumberOfParticles) {
  int numPtclsOld = blRingBufferSize(ensemble.buffer);
  blEnsembleRemoveBelow(0.5, ensemble.x, &ensemble);
  int numPtclsNew = blRingBufferSize(ensemble.buffer);
  assert_that(numPtclsOld >= numPtclsNew, is_true);
}

Ensure(RemoveBelow, removesParticlesBelow) {
  int i;
  blEnsembleRemoveBelow(0.5, ensemble.x, &ensemble);
  for (i = ensemble.buffer.begin; i != ensemble.buffer.end;
       i = blRingBufferNext(ensemble.buffer, i)) {
    assert_that_double(ensemble.x[i], is_greater_than_double(0.4999));
  }
}

Ensure(RemoveBelow, leavesParticlesAbove) {
  struct BLRingBuffer buffOld = ensemble.buffer;
  blEnsembleRemoveBelow(0.5, ensemble.x, &ensemble);
  int i;
  for (i = buffOld.begin; i != ensemble.buffer.begin;
       i = blRingBufferNext(buffOld, i)) {
    assert_that_double(ensemble.x[i], is_less_than_double(0.500001));
  }
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, RemoveBelow, doesNotIncreaseNumberOfParticles);
  add_test_with_context(suite, RemoveBelow, removesParticlesBelow);
  add_test_with_context(suite, RemoveBelow, leavesParticlesAbove);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

