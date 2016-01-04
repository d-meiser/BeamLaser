#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <Errors.h>
#include <Utilities.h>


struct BLEnsemble {
  int numPtcls;
  int maxNumPtcls;
  int internalStateSize;
  double *x;
  double *y;
  double *z;
  double *vx;
  double *vy;
  double *vz;
  double *internalState;
};

BL_STATUS blEnsembleInitialize(int capacity, int internalStateSize,
                               struct BLEnsemble *ensemble);
void blEnsembleFree(struct BLEnsemble *ensemble);
void blEnsembleRemoveBelow(double cutoff, double *positions,
                           struct BLEnsemble *ensemble);
void blEnsemblePush(double dt, struct BLEnsemble *ensemble);

#endif

