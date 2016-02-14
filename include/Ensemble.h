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

BL_STATUS blEnsembleCreate(int capacity, int internalStateSize,
                               struct BLEnsemble *ensemble);
void blEnsembleDestroy(struct BLEnsemble *ensemble);
void blEnsembleRemoveBelow(double cutoff, double *positions,
                           struct BLEnsemble *ensemble);
void blEnsembleCreateSpace(int numParticles, struct BLEnsemble *ensemble);
void blEnsemblePush(double dt, struct BLEnsemble *ensemble);

#endif

