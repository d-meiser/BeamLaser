#include <BeamLaser.h>
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define H_BAR 1.0545718e-34
#define EPSILON_0 8.85e-12
#define SPEED_OF_LIGHT 299792458.0

int main() {
  int maxNumPtcls = 100;
  int internalStateSize = 2;
  double complex initialState[2];
  initialState[0] = 0;
  initialState[1] = 1;

  struct BLSimulationState simulationState;
  blEnsembleCreate(maxNumPtcls, internalStateSize,
      &simulationState.ensemble);
  simulationState.ensemble.ptclWeight = 0;
  simulationState.fieldState.q = 1.0;
  simulationState.fieldState.p = 0.0;

  double lambda = 1.0e-6;
  struct ParticleSource* src = blParticleSourceManualCreate(0.25*lambda, 0, 0, 0, 0, 0,
      internalStateSize, initialState, 0);
  blEnsembleCreateSpace(1, &simulationState.ensemble);
  blParticleSourceCreateParticles(src,
      simulationState.ensemble.x,
      simulationState.ensemble.y,
      simulationState.ensemble.z,
      simulationState.ensemble.vx,
      simulationState.ensemble.vy,
      simulationState.ensemble.vz,
      internalStateSize,
      simulationState.ensemble.internalState);

  double dipoleMatrixElement = 1.0e-29;
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
  printf("omega == %e\n", omega);

  struct BLAtomFieldInteraction *atomFieldInteraction =
    blAtomFieldInteractionCreate(
      simulationState.ensemble.maxNumPtcls,
      simulationState.ensemble.internalStateSize, 
      dipoleOperator, modeFunction);

  double dt = 1.0e-2 / omega;
  int i;
  int numSteps = 100;
  int dumpPeriod = 10;
  for (i = 0; i < numSteps; ++i) {
    if (i % dumpPeriod == 0) {
      printf("%4d  %e %e  %e %e  %e\n", i,
          creal(simulationState.ensemble.internalState[0]),
          cimag(simulationState.ensemble.internalState[0]),
          creal(simulationState.ensemble.internalState[1]),
          cimag(simulationState.ensemble.internalState[1]),
          cos(omega * i * dt));
    }
    blAtomFieldInteractionTakeStep(atomFieldInteraction,
        dt, &simulationState.fieldState, &simulationState.ensemble);
  }
  assert(fabs(
        creal(simulationState.ensemble.internalState[1]) -
        cos(omega * dt * numSteps)) < 1.0e-6);

  blAtomFieldInteractionDestroy(atomFieldInteraction);
  blModeFunctionDestroy(modeFunction);
  blDipoleOperatorDestroy(dipoleOperator);
  blEnsembleDestroy(&simulationState.ensemble);
  blParticleSourceDestroy(src);
  return 0;
}

