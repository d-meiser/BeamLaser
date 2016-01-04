#include <cgreen/cgreen.h>
#include <ParticleSource.h>

#define INTERNAL_STATE_SIZE 4
#define ENSEMBLE_CAPACITY 10
static double x[ENSEMBLE_CAPACITY];
static double y[ENSEMBLE_CAPACITY];
static double z[ENSEMBLE_CAPACITY];
static double vx[ENSEMBLE_CAPACITY];
static double vy[ENSEMBLE_CAPACITY];
static double vz[ENSEMBLE_CAPACITY];
static double internalState[ENSEMBLE_CAPACITY * INTERNAL_STATE_SIZE];

Describe(ParticleSource)
BeforeEach(ParticleSource) {}
AfterEach(ParticleSource) {}

Ensure(ParticleSource, nullSourceCreatesZeroParticles) {
  int numPtcls = blParticleSourceGetNumParticles(0);
  assert_that(numPtcls, is_equal_to(0));
}

Ensure(ParticleSource, uniformSourceCreatesParticlesInBox) {
  struct BBox volume = {0.0, 1.0, -1.0, 0.5, -2.5, 3.0};
  int numPtcls = 5;
  double vbar[3] = {0.0, 0.0, -200.0};
  double deltaV[3] = {1.0, 1.0, 10.0};
  double initialState[4] = {0, 0, 1.0, 0};
  int i;

  struct ParticleSource *particleSource = blParticleSourceUniformCreate(
      volume, numPtcls, vbar, deltaV, INTERNAL_STATE_SIZE, initialState, 0);

  blParticleSourceCreateParticles(particleSource, x, y, z, vx, vy, vz,
      internalState);

  for (i = 0; i < numPtcls; ++i) {
    assert_that_double(x[i], is_greater_than_double(volume.xmin - 1.0e-6));
    assert_that_double(x[i], is_less_than_double(volume.xmax + 1.0e-6));
    assert_that_double(y[i], is_greater_than_double(volume.ymin - 1.0e-6));
    assert_that_double(y[i], is_less_than_double(volume.ymax + 1.0e-6));
    assert_that_double(z[i], is_greater_than_double(volume.zmin - 1.0e-6));
    assert_that_double(z[i], is_less_than_double(volume.zmax + 1.0e-6));
  }
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, ParticleSource, nullSourceCreatesZeroParticles);
  add_test_with_context(suite, ParticleSource, uniformSourceCreatesParticlesInBox);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

