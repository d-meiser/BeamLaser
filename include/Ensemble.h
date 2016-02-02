#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <Errors.h>
#include <Utilities.h>
#include <complex.h>


struct BLEnsemble {
  int numPtcls;
  int maxNumPtcls;
  int internalStateSize;
  double ptclWeight;
  double *x;
  double *y;
  double *z;
  double *vx;
  double *vy;
  double *vz;
  double complex *internalState;
};

BL_STATUS blEnsembleInitialize(int capacity, int internalStateSize,
                               struct BLEnsemble *ensemble);
void blEnsembleFree(struct BLEnsemble *ensemble);
void blEnsembleRemoveBelow(double cutoff, double *positions,
                           struct BLEnsemble *ensemble);
void blEnsembleCreateSpace(int numParticles, struct BLEnsemble *ensemble);
void blEnsemblePush(double dt, struct BLEnsemble *ensemble);

#endif

