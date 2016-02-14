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

static struct BLParticleSource *particleSource;
static double complex internalState[2];

Describe(ManualSource)
BeforeEach(ManualSource) {
  internalState[0] = 0;
  internalState[1] = 1;
}
AfterEach(ManualSource) {}


Ensure(ManualSource, canBeCreated) {
  particleSource = blParticleSourceManualCreate(0, 0, 0, 0, 0, 0,
      2, internalState, 0);
  assert_that(particleSource, is_not_null);
  blParticleSourceDestroy(particleSource);
}

Ensure(ManualSource, createsOneParticle) {
  particleSource = blParticleSourceManualCreate(0, 0, 0, 0, 0, 0,
      2, internalState, 0);
  int n = blParticleSourceGetNumParticles(particleSource);
  assert_that(n, is_equal_to(1));
  blParticleSourceDestroy(particleSource);
}

Ensure(ManualSource, putsParticleSourceInCorrectLocation) {
  particleSource = blParticleSourceManualCreate(0.3, 0.2, -0.3, 1.4, 2.7, -1.0,
      2, internalState, 0);
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  blParticleSourceCreateParticles(particleSource,
      &x, &y, &z, &vx, &vy, &vz, 2, internalState);
  assert_that_double(x, is_equal_to_double(0.3));
  assert_that_double(y, is_equal_to_double(0.2));
  assert_that_double(z, is_equal_to_double(-0.3));
  assert_that_double(vx, is_equal_to_double(1.4));
  assert_that_double(vy, is_equal_to_double(2.7));
  assert_that_double(vz, is_equal_to_double(-1.0));
  assert_that_double(cabs(internalState[0]), is_equal_to_double(0.0));
  blParticleSourceDestroy(particleSource);
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, ManualSource, canBeCreated);
  add_test_with_context(suite, ManualSource, createsOneParticle);
  add_test_with_context(suite, ManualSource, putsParticleSourceInCorrectLocation);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

