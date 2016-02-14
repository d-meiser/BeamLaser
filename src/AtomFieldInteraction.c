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


struct BLAtomFieldInteraction* blAtomFieldInteractionCreate(int maxNumParticles,
    int internalStateSize,
    struct BLDipoleOperator *dipoleOperator,
    struct BLModeFunction *modeFunction) {
  struct BLAtomFieldInteraction *this = malloc(sizeof(*this));
  this->dipoleOperator = dipoleOperator;
  this->modeFunction = modeFunction;
  this->ex = malloc(maxNumParticles * sizeof(*this->ex));
  this->ey = malloc(maxNumParticles * sizeof(*this->ey));
  this->ez = malloc(maxNumParticles * sizeof(*this->ez));
  this->dx = malloc(maxNumParticles * sizeof(*this->dx));
  this->dy = malloc(maxNumParticles * sizeof(*this->dy));
  this->dz = malloc(maxNumParticles * sizeof(*this->dz));
  this->phix = malloc(maxNumParticles * sizeof(*this->phix));
  this->phiy = malloc(maxNumParticles * sizeof(*this->phiy));
  this->phiz = malloc(maxNumParticles * sizeof(*this->phiz));
  blIntegratorCreate("RK4", maxNumParticles *
                     sizeof(double complex) / sizeof(double) *
                     internalStateSize,
                     &this->integrator);
  return this;
}

void blAtomFieldInteractionDestroy(struct BLAtomFieldInteraction* atomFieldInteraction) {
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

void blAtomFieldInteractionTakeStep(struct BLAtomFieldInteraction *atomFieldInteraction,
                                    double dt,
                                    struct BLFieldState *fieldState,
                                    struct BLEnsemble *ensemble) {
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
