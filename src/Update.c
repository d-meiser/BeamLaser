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
#include <Update.h>
#include <stdlib.h>

void blUpdateTakeStep(struct BLUpdate *update, double t, double dt,
    struct BLSimulationState *state) {
  update->takeStep(t, dt, state, update->ctx);
}

void blUpdateDestroy(struct BLUpdate *update) {
  update->destroy(update->ctx);
  free(update);
}

static void blUpdateIdentityTakeStep(
    double t, double dt, struct BLSimulationState *state, void *ctx) {
  BL_UNUSED(t);
  BL_UNUSED(dt);
  BL_UNUSED(state);
  BL_UNUSED(ctx);
}

static void blUpdateIdentityDestroy(void *ctx) {
  BL_UNUSED(ctx);
}

struct BLUpdate *blUpdateIdentityCreate() {
  struct BLUpdate *this = malloc(sizeof(*this));
  this->takeStep = blUpdateIdentityTakeStep;
  this->destroy = blUpdateIdentityDestroy;
  this->ctx = 0;
  return this;
}

