#include <cgreen/cgreen.h>
#include <Ensemble.h>

static struct BLEnsemble ensemble;
static struct BBox box = {-0.1, 0.1, -0.2, 0.2, 0.3, 0.4};
static const int numPtcls = 10;

Describe(CreateParticle)
BeforeEach(CreateParticle) {
  blEnsembleInitialize(numPtcls + 1, 2, &ensemble);
  ensemble.buffer.end = numPtcls;
}
AfterEach(CreateParticle) {
  blEnsembleFree(&ensemble);
}

Ensure(CreateParticle, positionsWithinBoundingBox) {
  int i;
  for (i = 0; i < numPtcls; ++i) {
    blEnsembleCreateParticle(box, 100.0, 10.0, 1.0, i, &ensemble);
  }
  for (i = 0; i < numPtcls; ++i) {
    assert_that_double(ensemble.x[i],
                       is_greater_than_double(box.xmin - 1.0e-6));
    assert_that_double(ensemble.x[i],
                       is_less_than_double(box.xmax + 1.0e-6));
    assert_that_double(ensemble.y[i],
                       is_greater_than_double(box.ymin - 1.0e-6));
    assert_that_double(ensemble.y[i],
                       is_less_than_double(box.ymax + 1.0e-6));
    assert_that_double(ensemble.z[i],
                       is_greater_than_double(box.zmin - 1.0e-6));
    assert_that_double(ensemble.z[i],
                       is_less_than_double(box.zmax + 1.0e-6));
  }
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, CreateParticle, positionsWithinBoundingBox);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

