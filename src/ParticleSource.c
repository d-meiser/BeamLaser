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
#include <ParticleSource.h>
#include <BeamLaserConfig.h>
#include <stdlib.h>
#include <string.h>
#define SIMPLE_SPRNG
#include <sprng.h>


int blParticleSourceGetNumParticles(struct BLParticleSource *particleSource) {
  if (particleSource) {
    return particleSource->getNumParticles(particleSource->ctx) +
      blParticleSourceGetNumParticles(particleSource->next);
  } else {
    return 0;
  }
}

void blParticleSourceCreateParticles(struct BLParticleSource *particleSource,
    double *x, double *y, double *z, double *vx, double *vy, double *vz,
    int internalStateSize, double complex *internalState) {
  if (particleSource) {
    particleSource->createParticles(x, y, z, vx, vy, vz, internalState,
        particleSource->ctx);
    int n = particleSource->getNumParticles(particleSource->ctx);
    blParticleSourceCreateParticles(particleSource->next,
        x + n, y + n, z + n, vx + n, vy + n, vz + n,
        internalStateSize, internalState + n * internalStateSize);
  }
}

/* Create a ParticleSource node; concrete implementation has to finish
 * initialization */
static void blParticleSourceCreate(struct BLParticleSource *next,
                                   struct BLParticleSource **particleSource) {
  *particleSource = malloc(sizeof(**particleSource));
  (*particleSource)->next = next;
}

void blParticleSourceDestroy(struct BLParticleSource *particleSource) {
  if (particleSource) {
    struct BLParticleSource *next = particleSource->next;
    particleSource->destroy(particleSource->ctx);
    free(particleSource);
    blParticleSourceDestroy(next);
  }
}

/*
 * Implementation of the uniform particle source
 */
struct UniformCtx {
  struct BlBox volume;
  double nbar;
  int nextNumPtcls;
  double vbar[3];
  double deltaV[3];
  int internalStateSize;
  double complex *internalState;
};

static int uniformGetNumParticles(void *ctx) {
  struct UniformCtx *uctx = ctx;
  if (uctx->nextNumPtcls < 0) {
    uctx->nextNumPtcls = blGeneratePoisson(uctx->nbar);
  }
  return uctx->nextNumPtcls;
}

static void uniformCreateParticles(double *x, double *y, double *z,
      double *vx, double *vy, double *vz,
      double complex *internalState, void *ctx) {
  struct UniformCtx *uctx = ctx;
  const struct BlBox *box = &uctx->volume;
  int i;

  if (uctx->nextNumPtcls < 0) {
    uctx->nextNumPtcls = blGeneratePoisson(uctx->nbar);
  }

  for (i = 0; i < uctx->nextNumPtcls; ++i) {
    x[i] = box->xmin + (box->xmax - box->xmin) * sprng();
    y[i] = box->ymin + (box->ymax - box->ymin) * sprng();
    z[i] = box->zmin + (box->zmax - box->zmin) * sprng();
    vx[i] = blGenerateGaussianNoise(uctx->vbar[0], uctx->deltaV[0]);
    vy[i] = blGenerateGaussianNoise(uctx->vbar[1], uctx->deltaV[1]);
    vz[i] = blGenerateGaussianNoise(uctx->vbar[2], uctx->deltaV[2]);
    memcpy(&internalState[i * uctx->internalStateSize],
        uctx->internalState,
        uctx->internalStateSize * sizeof(*uctx->internalState));
  }

  uctx->nextNumPtcls = -1;
}

void uniformDestroy(void *ctx) {
  struct UniformCtx *uctx = ctx;
  free(uctx->internalState);
  free(uctx);
}

struct BLParticleSource *blParticleSourceUniformCreate(
    struct BlBox volume, double nbar, double *vbar, double *deltaV,
    int internalStateSize, double complex *internalState,
    struct BLParticleSource *next) {
  struct BLParticleSource *self;
  blParticleSourceCreate(next, &self);
  self->getNumParticles = uniformGetNumParticles;
  self->createParticles = uniformCreateParticles;
  self->destroy = uniformDestroy;
  struct UniformCtx *ctx = malloc(sizeof(struct UniformCtx));
  ctx->volume = volume;
  ctx->nbar = nbar;
  ctx->nextNumPtcls = -1;
  ctx->vbar[0] = vbar[0];
  ctx->vbar[1] = vbar[1];
  ctx->vbar[2] = vbar[2];
  ctx->deltaV[0] = deltaV[0];
  ctx->deltaV[1] = deltaV[1];
  ctx->deltaV[2] = deltaV[2];
  ctx->internalStateSize = internalStateSize;
  ctx->internalState = malloc(internalStateSize * sizeof(*ctx->internalState));
  memcpy(ctx->internalState, internalState,
      internalStateSize * sizeof(*ctx->internalState));
  self->ctx = ctx;
  return self;
}


/* Manual particle source */

struct ManualCtx {
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  int internalStateSize;
  double complex *internalState;
};

static int manualGetNumParticles(void *ctx) {
  BL_UNUSED(ctx);
  return 1;
}

static void manualCreateParticles(double *x, double *y, double *z,
      double *vx, double *vy, double *vz,
      double complex *internalState, void *c) {
  struct ManualCtx *ctx = (struct ManualCtx*)c;
  x[0] = ctx->x;
  y[0] = ctx->y;
  z[0] = ctx->z;
  vx[0] = ctx->vx;
  vy[0] = ctx->vy;
  vz[0] = ctx->vz;
  memcpy(internalState, ctx->internalState, ctx->internalStateSize * sizeof(*ctx->internalState));
}

void manualDestroy(void *c) {
  struct ManualCtx *ctx = c;
  free(ctx->internalState);
  free(ctx);
}

struct BLParticleSource *blParticleSourceManualCreate(
    double x, double y, double z, double vx, double vy, double vz,
    int internalStateSize, double complex *internalState,
    struct BLParticleSource *next) {
  struct BLParticleSource *self;
  blParticleSourceCreate(next, &self);
  self->getNumParticles = manualGetNumParticles;
  self->createParticles = manualCreateParticles;
  self->destroy = manualDestroy;
  struct ManualCtx *ctx = malloc(sizeof(*ctx));
  ctx->x = x;
  ctx->y = y;
  ctx->z = z;
  ctx->vx = vx;
  ctx->vy = vy;
  ctx->vz = vz;
  ctx->internalStateSize = internalStateSize;
  ctx->internalState = malloc(internalStateSize * sizeof(*ctx->internalState));
  memcpy(ctx->internalState, internalState,
      internalStateSize * sizeof(*ctx->internalState));
  self->ctx = ctx;
  return self;
}


