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
#include <AtomFieldInteraction.h>
#include <stdlib.h>
#include <Integrator.h>
#include <SimulationState.h>
#include <Update.h>


struct BLAtomFieldInteraction {
  struct BLDipoleOperator *dipoleOperator;
  struct BLModeFunction *modeFunction;
  double complex *ex;
  double complex *ey;
  double complex *ez;
  double complex *dx;
  double complex *dy;
  double complex *dz;
  double complex *phix;
  double complex *phiy;
  double complex *phiz;
  BLIntegrator integrator;
};

struct AtomFieldInteractionIntegratorContext {
  struct BLAtomFieldInteraction *atomFieldInteraction;
  struct BLEnsemble *ensemble;
};

static void interactionRHS(double t, int n, const double *x, double *y,
                           void *ctx);


static void blAtomFieldInteractionDestroy(
    void *ctx) {
  struct BLAtomFieldInteraction *atomFieldInteraction =
    (struct BLAtomFieldInteraction*)ctx;
  free(atomFieldInteraction->ex);
  free(atomFieldInteraction->ey);
  free(atomFieldInteraction->ez);
  free(atomFieldInteraction->dx);
  free(atomFieldInteraction->dy);
  free(atomFieldInteraction->dz);
  free(atomFieldInteraction->phix);
  free(atomFieldInteraction->phiy);
  free(atomFieldInteraction->phiz);
  blIntegratorDestroy(&atomFieldInteraction->integrator);
  free(atomFieldInteraction);
}

static void blAtomFieldInteractionTakeStep(
    double t,
    double dt,
    struct BLSimulationState *state,
    void *ctx) {
  BL_UNUSED(t);
  struct BLAtomFieldInteraction *atomFieldInteraction = 
    (struct BLAtomFieldInteraction*)ctx;
  struct BLEnsemble *ensemble = &state->ensemble;
  struct BLFieldState *fieldState = &state->fieldState;
  struct BLModeFunction *modeFunction = atomFieldInteraction->modeFunction;
  const int numPtcls = ensemble->numPtcls;
  const int fieldOffset = numPtcls * ensemble->internalStateSize;
  int n = 2 * numPtcls * ensemble->internalStateSize + 2;

  /* Pack field into internal state buffer. Note that the internalState array
  needs to have room for at least two additional doubles. After the field has
  been scattered to each rank we integrate its equations of motion redundantly.
  */
  BL_MPI_Request fieldRequest =
    blBcastBegin((const double*)fieldState,
        (double*)(ensemble->internalState + fieldOffset), 2);
  blModeFunctionEvaluate(modeFunction, ensemble->numPtcls,
      ensemble->x, ensemble->y, ensemble->z,
      atomFieldInteraction->phix,
      atomFieldInteraction->phiy,
      atomFieldInteraction->phiz);
  blBcastEnd(fieldRequest,
      (const double*)fieldState,
      (double*)(ensemble->internalState + fieldOffset), 2);

  struct AtomFieldInteractionIntegratorContext integratorCtx;
  integratorCtx.atomFieldInteraction = atomFieldInteraction;
  integratorCtx.ensemble = ensemble;
  blIntegratorTakeStep(atomFieldInteraction->integrator,
      0.0, dt, n, interactionRHS,
      (const double*)ensemble->internalState,
      (double*)ensemble->internalState,
      &integratorCtx);

  fieldState->q = creal(ensemble->internalState[fieldOffset]);
  fieldState->p = cimag(ensemble->internalState[fieldOffset]);
}

void interactionRHS(double t, int n, const double *x, double *y,
                    void *c) {
  BL_UNUSED(t);
  BL_UNUSED(n);
  struct AtomFieldInteractionIntegratorContext *ctx =
    (struct AtomFieldInteractionIntegratorContext*)c;
  struct BLAtomFieldInteraction *atomFieldInteraction = ctx->atomFieldInteraction;
  struct BLEnsemble *ensemble = ctx->ensemble;
  int i;
  const int numPtcls = ensemble->numPtcls;
  const int fieldOffset = numPtcls * ensemble->internalStateSize;

  blDipoleOperatorComputeD(atomFieldInteraction->dipoleOperator,
      ensemble->internalStateSize, numPtcls, (const double complex*)x,
      atomFieldInteraction->dx,
      atomFieldInteraction->dy,
      atomFieldInteraction->dz);

  double complex polarization = 0;
  for (i = 0; i < numPtcls; ++i) {
    polarization += conj(atomFieldInteraction->phix[i]) * atomFieldInteraction->dx[i];
  }
  for (i = 0; i < numPtcls; ++i) {
    polarization += conj(atomFieldInteraction->phiy[i]) * atomFieldInteraction->dy[i];
  }
  for (i = 0; i < numPtcls; ++i) {
    polarization += conj(atomFieldInteraction->phiz[i]) * atomFieldInteraction->dz[i];
  }
  polarization *= -I * ensemble->ptclWeight;

  BL_MPI_Request polReq =
    blAddAllBegin((const double*)&polarization, y + 2 * fieldOffset, 2);

  const double complex fieldAmplitude =
    *((const double complex*)(x + 2 * fieldOffset));
  for (i = 0; i < numPtcls; ++i) {
    atomFieldInteraction->ex[i] = fieldAmplitude * atomFieldInteraction->phix[i];
  }
  for (i = 0; i < numPtcls; ++i) {
    atomFieldInteraction->ey[i] = fieldAmplitude * atomFieldInteraction->phiy[i];
  }
  for (i = 0; i < numPtcls; ++i) {
    atomFieldInteraction->ez[i] = fieldAmplitude * atomFieldInteraction->phiz[i];
  }

  double complex *yc = (double complex*)y;
  blDipoleOperatorApply(atomFieldInteraction->dipoleOperator,
                        ensemble->internalStateSize,
                        numPtcls,
                        atomFieldInteraction->ex, atomFieldInteraction->ey, atomFieldInteraction->ez,
                        (const double complex*)x, yc);
  for (i = 0; i < numPtcls * ensemble->internalStateSize; ++i) {
    yc[i] *= -I;
  }

  blAddAllEnd(polReq, (const double*)&polarization, y + 2 * fieldOffset, 2);
}

struct BLUpdate* blAtomFieldInteractionCreate(int maxNumParticles,
    int internalStateSize, struct BLDipoleOperator *dipoleOperator,
    struct BLModeFunction *modeFunction) {
  struct BLUpdate *this = malloc(sizeof(*this));
  this->takeStep = blAtomFieldInteractionTakeStep;
  this->destroy = blAtomFieldInteractionDestroy;
  struct BLAtomFieldInteraction *ctx = malloc(sizeof(*ctx));
  ctx->dipoleOperator = dipoleOperator;
  ctx->modeFunction = modeFunction;
  ctx->ex = malloc(maxNumParticles * sizeof(*ctx->ex));
  ctx->ey = malloc(maxNumParticles * sizeof(*ctx->ey));
  ctx->ez = malloc(maxNumParticles * sizeof(*ctx->ez));
  ctx->dx = malloc(maxNumParticles * sizeof(*ctx->dx));
  ctx->dy = malloc(maxNumParticles * sizeof(*ctx->dy));
  ctx->dz = malloc(maxNumParticles * sizeof(*ctx->dz));
  ctx->phix = malloc(maxNumParticles * sizeof(*ctx->phix));
  ctx->phiy = malloc(maxNumParticles * sizeof(*ctx->phiy));
  ctx->phiz = malloc(maxNumParticles * sizeof(*ctx->phiz));
  blIntegratorCreate("RK4", maxNumParticles *
                     sizeof(double complex) / sizeof(double) *
                     internalStateSize,
                     &ctx->integrator);
  this->ctx = ctx;
  return this;
}

