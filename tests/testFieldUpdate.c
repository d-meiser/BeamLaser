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
#include <FieldUpdate.h>
#include <math.h>

static struct BLUpdate *update;

Describe(FieldUpdate)
BeforeEach(FieldUpdate) {}
AfterEach(FieldUpdate) {}


Ensure(FieldUpdate, canBeCreated) {
  update = blFieldUpdateCreate(0.0, 1.0, 0.0);
  assert_that(update, is_not_null);
  blUpdateDestroy(update);
}

Ensure(FieldUpdate, givesExponentialDamping) {
  double kappa = 3.8;
  update = blFieldUpdateCreate(0.0, kappa, 0.0);
  struct BLSimulationState state;
  state.fieldState.q = 2.3;
  state.fieldState.p = -1.8;
  double dt = 1.7;
  blUpdateTakeStep(update, 0.0, dt, &state);
  assert_that_double(state.fieldState.q, is_equal_to_double(2.3 * exp(-kappa * dt)));
  assert_that_double(state.fieldState.p, is_equal_to_double(-1.8 * exp(-kappa * dt)));
  blUpdateDestroy(update);
}

Ensure(FieldUpdate, givesPhaseRotation) {
  double delta = 3.8;
  update = blFieldUpdateCreate(delta, 0.0, 0.0);
  struct BLSimulationState state;
  state.fieldState.q = 2.3;
  state.fieldState.p = -1.8;
  double dt = 1.7;
  blUpdateTakeStep(update, 0.0, dt, &state);
  assert_that_double(state.fieldState.q,
      is_equal_to_double(2.3 * cos(delta * dt) + 1.8 * sin(delta * dt)));
  assert_that_double(state.fieldState.p,
    is_equal_to_double(2.3 * sin(delta * dt) - 1.8 * cos(delta * dt)));
  blUpdateDestroy(update);
}

Ensure(FieldUpdate, givesBoundedEvolutionWithNoise) {
  double kappa = 3.8;
  update = blFieldUpdateCreate(0.0, kappa, sqrt(kappa));
  struct BLSimulationState state;
  state.fieldState.q = 2.3;
  state.fieldState.p = -1.8;
  double dt = 1.7;
  int i;
  for (i = 0; i < 10000; ++i) {
    blUpdateTakeStep(update, i * dt, dt, &state);
  }
  assert_that_double(fabs(state.fieldState.q), is_less_than_double(10.0));
  assert_that_double(fabs(state.fieldState.p), is_less_than_double(10.0));
  blUpdateDestroy(update);
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
  add_test_with_context(suite, FieldUpdate, canBeCreated);
  add_test_with_context(suite, FieldUpdate, givesExponentialDamping);
  add_test_with_context(suite, FieldUpdate, givesPhaseRotation);
  add_test_with_context(suite, FieldUpdate, givesBoundedEvolutionWithNoise);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return result;
}


