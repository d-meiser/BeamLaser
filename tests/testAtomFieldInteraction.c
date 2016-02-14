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
#include <AtomFieldInteraction.h>
#include <SimulationState.h>
#include <Update.h>
#include <math.h>


#define MAX_NUM_PTCLS 10
#define DOF_PER_PTCL 2

static struct BLDipoleOperator *dipoleOperator;
static struct BLModeFunction *modeFunction;
static struct BLSimulationState simulationState;
static struct BLUpdate *atomFieldInteraction;
static int i, j;

Describe(AtomFieldInteraction)

BeforeEach(AtomFieldInteraction) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  modeFunction = blModeFunctionUniformCreate(0.0, 1.0, 0.0);

  simulationState.fieldState.q = 1.0;
  simulationState.fieldState.p = 0.0;

  blEnsembleCreate(MAX_NUM_PTCLS, DOF_PER_PTCL, &simulationState.ensemble);
  simulationState.ensemble.numPtcls = 1;
  for (i = 0; i < MAX_NUM_PTCLS; ++i) {
    simulationState.ensemble.x[i] = 0;
    simulationState.ensemble.y[i] = 0;
    simulationState.ensemble.z[i] = 0;
    simulationState.ensemble.vx[i] = 0;
    simulationState.ensemble.vy[i] = 0;
    simulationState.ensemble.vz[i] = 0;
    for (j = 0; j < DOF_PER_PTCL; ++j) {
      simulationState.ensemble.internalState[i * DOF_PER_PTCL + j] = (j == 1) ? 1.0 : 0.0;
    }
  }
  atomFieldInteraction =
    blAtomFieldInteractionCreate(MAX_NUM_PTCLS, DOF_PER_PTCL, dipoleOperator, modeFunction);
}

AfterEach(AtomFieldInteraction) {
  blEnsembleDestroy(&simulationState.ensemble);
  blUpdateDestroy(atomFieldInteraction);
  blModeFunctionDestroy(modeFunction);
  blDipoleOperatorDestroy(dipoleOperator);
}


Ensure(AtomFieldInteraction, canBeCreated) {
  assert_that(atomFieldInteraction, is_not_null);
}

static double nrm_squared(double complex z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

Ensure(AtomFieldInteraction, isNormConserving) {
  for (i = 0; i < 1000; ++i) {
    blUpdateTakeStep(atomFieldInteraction,
        i * 1.0e-3, 1.0e-3, &simulationState);
  }
  double nrm = 0;
  for (i = 0; i < 2; ++i) {
    nrm += nrm_squared(simulationState.ensemble.internalState[i]);
  }
  assert_that_double(nrm, is_equal_to_double(1.0));
}

Ensure(AtomFieldInteraction, producesRabiOscillations) {
  /* Use zero weight particles to produce a constant field */
  simulationState.ensemble.ptclWeight = 0.0;
  for (i = 0; i < 10000; ++i) {
    blUpdateTakeStep(atomFieldInteraction,
        i * 1.0e-3, 1.0e-3, &simulationState);
  }
  assert_that_double(simulationState.fieldState.q, is_equal_to_double(1.0));
  assert_that_double(simulationState.fieldState.p, is_equal_to_double(0.0));

  /* Have to check phase factors here */
  assert_that_double(creal(simulationState.ensemble.internalState[0]),
      is_equal_to_double(0.0));
  assert_that_double(cimag(simulationState.ensemble.internalState[0]),
      is_equal_to_double(-sin(10.0)));
  assert_that_double(creal(simulationState.ensemble.internalState[1]),
      is_equal_to_double(cos(10.0)));
  assert_that_double(cimag(simulationState.ensemble.internalState[1]),
      is_equal_to_double(0.0));
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
  add_test_with_context(suite, AtomFieldInteraction, canBeCreated);
  add_test_with_context(suite, AtomFieldInteraction, isNormConserving);
  add_test_with_context(suite, AtomFieldInteraction, producesRabiOscillations);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return result;
}

