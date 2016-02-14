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
#include <Source.h>
#include <stdlib.h>


struct SourceCtx {
  struct BLParticleSource *src;
};

static void sourceTakeStep(double t, double dt, struct BLSimulationState *state, void *c) {
  BL_UNUSED(t);
  BL_UNUSED(dt);
  struct SourceCtx *ctx = (struct SourceCtx *)c;
  int numCreate = blParticleSourceGetNumParticles(ctx->src);
  struct BLEnsemble *ensemble = &state->ensemble;
  blEnsembleCreateSpace(numCreate, ensemble);
  blParticleSourceCreateParticles(ctx->src,
                                  ensemble->x, ensemble->y, ensemble->z,
                                  ensemble->vx, ensemble->vy, ensemble->vz,
                                  ensemble->internalStateSize,
                                  ensemble->internalState);
}

static void sourceDestroy(void *c) {
  free(c);
}

struct BLUpdate *blSourceCreate(struct BLParticleSource *src) {
  struct BLUpdate *this = malloc(sizeof(*this));
  this->takeStep = sourceTakeStep;
  this->destroy = sourceDestroy;
  struct SourceCtx *ctx = malloc(sizeof(*ctx));
  ctx->src = src;
  this->ctx = ctx;
  return this;
}

