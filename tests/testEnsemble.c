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
#include <Ensemble.h>

static struct BLEnsemble ensemble;

Describe(Ensemble)
BeforeEach(Ensemble) {}
AfterEach(Ensemble) {}

Ensure(Ensemble, initiallyHasZeroSize) {
  blEnsembleCreate(10, 4, &ensemble);
  assert_that(ensemble.numPtcls, is_equal_to(0));
  blEnsembleDestroy(&ensemble);
}

Ensure(Ensemble, hasRequestedCapacity) {
  const int requestedCapacity = 6;
  blEnsembleCreate(requestedCapacity, 4, &ensemble);
  assert_that(ensemble.maxNumPtcls, is_equal_to(requestedCapacity));
  blEnsembleDestroy(&ensemble);
}

Ensure(Ensemble, canBeDestroyd) {
  blEnsembleCreate(5, 4, &ensemble);
  blEnsembleDestroy(&ensemble);
}

Ensure(Ensemble, createSpace) {
  blEnsembleCreate(5, 4, &ensemble);
  blEnsembleCreateSpace(2, &ensemble);
  assert_that(ensemble.numPtcls, is_equal_to(2));
  blEnsembleDestroy(&ensemble);
}

Ensure(Ensemble, growsWhenOutOfSpace) {
  blEnsembleCreate(5, 4, &ensemble);
  blEnsembleCreateSpace(10, &ensemble);
  int i, j;
  for (i = 0; i < 10; ++i) {
    ensemble.x[i] = 0;
    ensemble.y[i] = 0;
    ensemble.z[i] = 0;
    ensemble.vx[i] = 0;
    ensemble.vy[i] = 0;
    ensemble.vz[i] = 0;
    for (j = 0; j < ensemble.internalStateSize; ++j) {
      ensemble.internalState[i * ensemble.internalStateSize + j] = 0;
    }
  }
  blEnsembleDestroy(&ensemble);
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Ensemble, initiallyHasZeroSize);
  add_test_with_context(suite, Ensemble, hasRequestedCapacity);
  add_test_with_context(suite, Ensemble, canBeDestroyd);
  add_test_with_context(suite, Ensemble, createSpace);
  add_test_with_context(suite, Ensemble, growsWhenOutOfSpace);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}
