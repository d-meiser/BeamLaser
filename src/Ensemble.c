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
#include <Ensemble.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <Partition.h>


BL_STATUS blEnsembleCreate(int capacity, int internalStateSize,
                               struct BLEnsemble *ensemble) {
  ensemble->numPtcls = 0;
  ensemble->maxNumPtcls = capacity;

  ensemble->internalStateSize = internalStateSize;
  ensemble->ptclWeight = 1.0;

  ensemble->x = malloc(capacity * sizeof(*ensemble->x));
  if (!ensemble->x) { return BL_OUT_OF_MEMORY; }
  ensemble->y = malloc(capacity * sizeof(*ensemble->y));
  if (!ensemble->y) { return BL_OUT_OF_MEMORY; }
  ensemble->z = malloc(capacity * sizeof(*ensemble->z));
  if (!ensemble->z) { return BL_OUT_OF_MEMORY; }
  ensemble->vx = malloc(capacity * sizeof(*ensemble->vx));
  if (!ensemble->vx) { return BL_OUT_OF_MEMORY; }
  ensemble->vy = malloc(capacity * sizeof(*ensemble->vy));
  if (!ensemble->vy) { return BL_OUT_OF_MEMORY; }
  ensemble->vz = malloc(capacity * sizeof(*ensemble->vz));
  if (!ensemble->vz) { return BL_OUT_OF_MEMORY; }
  ensemble->internalState = malloc(capacity * internalStateSize *
                                   sizeof(*ensemble->internalState));
  if (!ensemble->internalState) { return BL_OUT_OF_MEMORY; }
  return BL_SUCCESS;
}

void blEnsembleDestroy(struct BLEnsemble *ensemble) {
  ensemble->numPtcls = 0;
  ensemble->maxNumPtcls = 0;

  free(ensemble->x);
  free(ensemble->y);
  free(ensemble->z);
  free(ensemble->vx);
  free(ensemble->vy);
  free(ensemble->vz);
  free(ensemble->internalState);
}

#define SWAP(a,b) do { \
  tmp = (a); \
  (a) = (b); \
  (b) = tmp; \
} while (0)

static void ensembleSwap_(int i, int j, void *ctx) {
  int m;
  struct BLEnsemble *ensemble = ctx;
  if (i != j) {
    double tmp;
    SWAP(ensemble->x[i], ensemble->x[j]);
    SWAP(ensemble->y[i], ensemble->y[j]);
    SWAP(ensemble->z[i], ensemble->z[j]);
    SWAP(ensemble->vx[i], ensemble->vx[j]);
    SWAP(ensemble->vy[i], ensemble->vy[j]);
    SWAP(ensemble->vz[i], ensemble->vz[j]);
    for (m = 0; m < ensemble->internalStateSize; ++m) {
      double complex ctmp =
        ensemble->internalState[ensemble->internalStateSize * i + m];
      ensemble->internalState[ensemble->internalStateSize * i + m] =
        ensemble->internalState[ensemble->internalStateSize * j + m];
      ensemble->internalState[ensemble->internalStateSize * j + m] = ctmp;
    }
  }
}

void blEnsembleRemoveBelow(double cutoff, double *positions,
                           struct BLEnsemble *ensemble) {
  struct BLSwapClosure swap;
  swap.f = ensembleSwap_;
  swap.ctx = ensemble;
  ensemble->numPtcls = blBSP(0, ensemble->numPtcls, cutoff, positions, swap);
}

void blEnsemblePush(double dt, struct BLEnsemble *ensemble) {
  int i;
  double * restrict x = ensemble->x;
  double * restrict y = ensemble->y;
  double * restrict z = ensemble->z;
  const double * restrict vx = ensemble->vx;
  const double * restrict vy = ensemble->vy;
  const double * restrict vz = ensemble->vz;
  for (i = 0; i < ensemble->numPtcls; ++i) {
    x[i] += dt * vx[i];
  }
  for (i = 0; i < ensemble->numPtcls; ++i) {
    y[i] += dt * vy[i];
  }
  for (i = 0; i < ensemble->numPtcls; ++i) {
    z[i] += dt * vz[i];
  }
}

void blEnsembleCreateSpace(int numParticles, struct BLEnsemble *ensemble) {
  memmove(ensemble->x + numParticles, ensemble->x,
      ensemble->numPtcls * sizeof(*ensemble->x));
  memmove(ensemble->y + numParticles, ensemble->y,
      ensemble->numPtcls * sizeof(*ensemble->y));
  memmove(ensemble->z + numParticles, ensemble->z,
      ensemble->numPtcls * sizeof(*ensemble->z));
  memmove(ensemble->vx + numParticles, ensemble->vx,
      ensemble->numPtcls * sizeof(*ensemble->vx));
  memmove(ensemble->vy + numParticles, ensemble->vy,
      ensemble->numPtcls * sizeof(*ensemble->vy));
  memmove(ensemble->vz + numParticles, ensemble->vz,
      ensemble->numPtcls * sizeof(*ensemble->vz));
  memmove(ensemble->internalState + numParticles * ensemble->internalStateSize,
      ensemble->internalState,
      ensemble->numPtcls * sizeof(*ensemble->internalState) *
      ensemble->internalStateSize);
  ensemble->numPtcls += numParticles;
}
