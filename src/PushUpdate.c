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
#include <PushUpdate.h>

static void pushUpdateTakeStep(double t, double dt, struct BLSimulationState *state, void *ctx) {
  BL_UNUSED(t);
  BL_UNUSED(ctx);
  int i;
  struct BLEnsemble *ensemble = &state->ensemble;
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

static void pushUpdateDestroy(void *ctx) {
  BL_UNUSED(ctx);
}

struct BLUpdate *blPushUpdateCreate() {
  struct BLUpdate *this = malloc(sizeof(*this));
  this->takeStep = pushUpdateTakeStep;
  this->destroy = pushUpdateDestroy;
  this->ctx = 0;
  return this;
}

