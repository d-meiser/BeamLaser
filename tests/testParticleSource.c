/*
Copyright 2014 Dominic Meiser

This file is part of BeamLaser.

BeamLaser is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

BeamLaser is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with BeamLaser.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <cgreen/cgreen.h>
#include <ParticleSource.h>
#include <config.h>
#ifdef BL_WITH_MPI
#include <mpi.h>
#endif

#define INTERNAL_STATE_SIZE 2
#define ENSEMBLE_CAPACITY 10
static double x[ENSEMBLE_CAPACITY];
static double y[ENSEMBLE_CAPACITY];
static double z[ENSEMBLE_CAPACITY];
static double vx[ENSEMBLE_CAPACITY];
static double vy[ENSEMBLE_CAPACITY];
static double vz[ENSEMBLE_CAPACITY];
static double complex internalState[ENSEMBLE_CAPACITY * INTERNAL_STATE_SIZE];
static struct BlBox volume = {0.0, 1.0, -1.0, 0.5, -2.5, 3.0};
static double vbar[3] = {0.0, 0.0, -200.0};
static double deltaV[3] = {1.0, 1.0, 10.0};
static double complex initialState[2] = {0, 1.0};
static const int numPtcls = 5;
static struct BLParticleSource *particleSource;

Describe(ParticleSource)
BeforeEach(ParticleSource) {
  particleSource = blParticleSourceUniformCreate(
      volume, numPtcls, vbar, deltaV, INTERNAL_STATE_SIZE, initialState, 0);
}

AfterEach(ParticleSource) {
  blParticleSourceDestroy(particleSource);
}

Ensure(ParticleSource, nullSourceCreatesZeroParticles) {
  int n = blParticleSourceGetNumParticles(0);
  assert_that(n, is_equal_to(0));
}

Ensure(ParticleSource, twoSourcesCanBeDestroyed) {
  particleSource = blParticleSourceUniformCreate(
      volume, numPtcls, vbar, deltaV, INTERNAL_STATE_SIZE, initialState,
      particleSource);
}

Ensure(ParticleSource, uniformSourceProducesNonNegativeNumberOfParticles) {
  int n = blParticleSourceGetNumParticles(particleSource);
  assert_that(n, is_greater_than(-1));
}

Ensure(ParticleSource, uniformSourceCreatesParticlesInBox) {
  int i;

  blParticleSourceCreateParticles(particleSource, x, y, z, vx, vy, vz,
      INTERNAL_STATE_SIZE, internalState);

  for (i = 0; i < numPtcls; ++i) {
    assert_that_double(x[i], is_greater_than_double(volume.xmin - 1.0e-6));
    assert_that_double(x[i], is_less_than_double(volume.xmax + 1.0e-6));
    assert_that_double(y[i], is_greater_than_double(volume.ymin - 1.0e-6));
    assert_that_double(y[i], is_less_than_double(volume.ymax + 1.0e-6));
    assert_that_double(z[i], is_greater_than_double(volume.zmin - 1.0e-6));
    assert_that_double(z[i], is_less_than_double(volume.zmax + 1.0e-6));
  }
}

Ensure(ParticleSource, twoSourcesCanBeComposed) {
  particleSource = blParticleSourceUniformCreate(
      volume, numPtcls, vbar, deltaV, INTERNAL_STATE_SIZE, initialState,
      particleSource);
  int i;
  for (i = 0; i < 2 * numPtcls; ++i) {
    assert_that_double(x[i], is_greater_than_double(volume.xmin - 1.0e-6));
    assert_that_double(x[i], is_less_than_double(volume.xmax + 1.0e-6));
    assert_that_double(y[i], is_greater_than_double(volume.ymin - 1.0e-6));
    assert_that_double(y[i], is_less_than_double(volume.ymax + 1.0e-6));
    assert_that_double(z[i], is_greater_than_double(volume.zmin - 1.0e-6));
    assert_that_double(z[i], is_less_than_double(volume.zmax + 1.0e-6));
  }
}

int main(int argn, char **argv)
{
#ifdef BL_WITH_MPI
  MPI_Init(&argn, &argv);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, ParticleSource, nullSourceCreatesZeroParticles);
  add_test_with_context(suite, ParticleSource, twoSourcesCanBeDestroyed);
  add_test_with_context(suite, ParticleSource,
      uniformSourceProducesNonNegativeNumberOfParticles);
  add_test_with_context(suite, ParticleSource, uniformSourceCreatesParticlesInBox);
  add_test_with_context(suite, ParticleSource, twoSourcesCanBeComposed);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return result;
}

