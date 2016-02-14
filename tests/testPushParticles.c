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
#include <PushUpdate.h>


static struct BLSimulationState simulationState;
static struct BLUpdate *push;

Describe(Push)
BeforeEach(Push) {
  push = blPushUpdateCreate();
  static const int numPtcls = 5;
  blEnsembleCreate(numPtcls + 1, 2, &simulationState.ensemble);
  int i;
  for (i = 0; i < numPtcls; ++i) {
    simulationState.ensemble.x[i] = 0;
    simulationState.ensemble.y[i] = 0;
    simulationState.ensemble.z[i] = 0;
    simulationState.ensemble.vx[i] = rand() / (double)RAND_MAX;
    simulationState.ensemble.vy[i] = rand() / (double)RAND_MAX;
    simulationState.ensemble.vz[i] = rand() / (double)RAND_MAX;
  }
  simulationState.ensemble.numPtcls = numPtcls;
}
AfterEach(Push) {
  blEnsembleDestroy(&simulationState.ensemble);
  blUpdateDestroy(push);
}

Ensure(Push, movesParticlesTheRightDistance) {
  blUpdateTakeStep(push, 0.0, 0.7, &simulationState);
  int i;
  for (i = 0; i != simulationState.ensemble.numPtcls; ++i) {
    assert_that_double(simulationState.ensemble.x[i],
        is_equal_to_double(simulationState.ensemble.vx[i] * 0.7));
    assert_that_double(simulationState.ensemble.y[i],
        is_equal_to_double(simulationState.ensemble.vy[i] * 0.7));
    assert_that_double(simulationState.ensemble.z[i],
        is_equal_to_double(simulationState.ensemble.vz[i] * 0.7));
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

