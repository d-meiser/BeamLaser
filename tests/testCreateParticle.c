#include <cgreen/cgreen.h>
#include <Ensemble.h>

static struct BLEnsemble ensemble;
static struct BBox box = {-0.1, 0.1, -0.2, 0.2, 0.3, 0.4};

Describe(CreateParticle)
BeforeEach(CreateParticle) {
  static const int numPtcls = 5;
  blEnsembleInitialize(numPtcls + 1, 2, &ensemble);
  ensemble.buffer.end = numPtcls;
}
AfterEach(CreateParticle) {
  blEnsembleFree(&ensemble);
}

Ensure(CreateParticle, positionsWithinBoundingBox) {
  blEnsembleCreateParticle(box, 100.0, 10.0, 1.0, 0, &ensemble);
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, CreateParticle, positionsWithinBoundingBox);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

