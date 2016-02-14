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
#include <BeamLaser.h>
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>


#ifdef BL_WITH_MPI
#include <mpi.h>
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define H_BAR 1.0545718e-34
#define EPSILON_0 8.85e-12
#define SPEED_OF_LIGHT 299792458.0


int main(int argn, char **argv) {
#ifdef BL_WITH_MPI
  MPI_Init(&argn, &argv);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif
  double nbar = 1.0e2;
  int maxNumPtcls = 100 * nbar;
  int internalStateSize = 2;
  double complex initialState[2];
  initialState[0] = 0;
  initialState[1] = 1;

  struct BLSimulationState simulationState;
  blEnsembleCreate(maxNumPtcls, internalStateSize,
      &simulationState.ensemble);
  simulationState.ensemble.ptclWeight = 100.0;
  simulationState.fieldState.q = 1.0;
  simulationState.fieldState.p = 0.0;

  double lambda = 1.0e-6;
  struct BlBox simulationBox =
    {-0.5 * lambda, 0.5 * lambda, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4};
  double vbar[3] = {0};
  double deltaV[3] = {0};
  struct ParticleSource* src =
    blParticleSourceUniformCreate(simulationBox, 1.0e4, vbar, deltaV,
                                  internalStateSize, initialState, 0);

  int n = blParticleSourceGetNumParticles(src);
  blEnsembleCreateSpace(n, &simulationState.ensemble);
  blParticleSourceCreateParticles(src,
      simulationState.ensemble.x,
      simulationState.ensemble.y,
      simulationState.ensemble.z,
      simulationState.ensemble.vx,
      simulationState.ensemble.vy,
      simulationState.ensemble.vz,
      internalStateSize,
      simulationState.ensemble.internalState);

  double dipoleMatrixElement = 1.0e-31;
  struct BLDipoleOperator *dipoleOperator =
    blDipoleOperatorTLACreate(dipoleMatrixElement);

  double waist = 1.0e-4;
  double length = 1.0e-2;
  struct BLModeFunction *modeFunction =
    blModeFunctionSimplifiedGaussianCreate(waist, lambda, length);

  double veff = waist * waist * length;
  double omega =
    sqrt(2.0 * M_PI * SPEED_OF_LIGHT /
        (lambda * 2.0 * EPSILON_0 * veff * H_BAR)) *
    dipoleMatrixElement;

  struct BLAtomFieldInteraction *atomFieldInteraction =
    blAtomFieldInteractionCreate(
      simulationState.ensemble.maxNumPtcls,
      simulationState.ensemble.internalStateSize, 
      dipoleOperator, modeFunction);

  double dt = 1.0e-3 / omega;
  printf("# dt == %le\n", dt);
  int i;
  int dumpPeriodicity = 10;
  struct BLDiagnostics *diagnostics = blDiagnosticsFieldStateCreate(dumpPeriodicity, 0);
  diagnostics = blDiagnosticsInternalStateCreate(dumpPeriodicity, "internal_state", diagnostics);
  int numSteps = 100;
  for (i = 0; i < numSteps; ++i) {
    blDiagnosticsProcess(diagnostics, i, &simulationState);
    blAtomFieldInteractionTakeStep(atomFieldInteraction,
        dt, &simulationState.fieldState, &simulationState.ensemble);
  }

  blDiagnosticsDestroy(diagnostics);
  blAtomFieldInteractionDestroy(atomFieldInteraction);
  blModeFunctionDestroy(modeFunction);
  blDipoleOperatorDestroy(dipoleOperator);
  blEnsembleDestroy(&simulationState.ensemble);
  blParticleSourceDestroy(src);

#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return BL_SUCCESS;
}

