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
#include <Sink.h>

struct SinkBelowCtx {
  double zmin;
};

static void sinkBelowTakeStep(double t, double dt, struct BLSimulationState *simulationState, void *c) {
  struct SinkBelowCtx *ctx = (struct SinkBelowCtx*)c;
  blEnsembleRemoveBelow(ctx->zmin, simulationState->ensemble.z, &simulationState->ensemble);
}

static void sinkBelowDestroy(void *c) {
  free(c);
}

struct BLUpdate *blSinkBelowCreate(double zmin) {
  struct BLUpdate *this = malloc(sizeof(*this));
  this->takeStep = sinkBelowTakeStep;
  this->destroy = sinkBelowDestroy;
  struct SinkBelowCtx *ctx = malloc(sizeof(*ctx));
  ctx->zmin = zmin;
  this->ctx = ctx;
  return this;
}

