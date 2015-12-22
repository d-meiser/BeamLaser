#include <cgreen/cgreen.h>
#include <Ensemble.h>

static struct BLEnsemble ensemble;

Describe(RemoveBelow)
BeforeEach(RemoveBelow) {
  blEnsembleInitialize(20, 2, &ensemble);
  int i;
  for (i = 0; i < 20; ++i) {
    ensemble.x[i] = rand() / (double)RAND_MAX;
  }
}
AfterEach(RemoveBelow) {
  blEnsembleFree(&ensemble);
}

Ensure(RemoveBelow, removesParticlesBelow) {
  blEnsembleRemoveBelow(0.5, ensemble.x, &ensemble);
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, RemoveBelow, removesParticlesBelow);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

