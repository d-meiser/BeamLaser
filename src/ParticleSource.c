#include <ParticleSource.h>
#include <stdlib.h>
#include <string.h>
#define SIMPLE_SPRNG
#include <sprng.h>

int blParticleSourceGetNumParticles(struct ParticleSource *particleSource) {
  if (particleSource) {
    return particleSource->getNumParticles(particleSource->ctx) +
      blParticleSourceGetNumParticles(particleSource->next);
  } else {
    return 0;
  }
}

void blParticleSourceCreateParticles(struct ParticleSource *particleSource,
    double *x, double *y, double *z, double *vx, double *vy, double *vz,
    double *internalState) {
  if (particleSource) {
    particleSource->createParticles(x, y, z, vx, vy, vz, internalState,
        particleSource->ctx);
    int n = particleSource->getNumParticles(particleSource->ctx);
    blParticleSourceCreateParticles(particleSource->next,
        x + n, y + n, z + n, vx + n, vy + n, vz + n, internalState + n);
  }
}

/* Create a ParticleSource node; concrete implementation has to finish
 * initialization */
static void blParticleSourceCreate(struct ParticleSource *next,
                                   struct ParticleSource **particleSource) {
  *particleSource = malloc(sizeof(**particleSource));
  (*particleSource)->next = next;
}

void blParticleSourceDestroy(struct ParticleSource *particleSource) {
  if (particleSource) {
    struct ParticleSource *next = particleSource->next;
    particleSource->destroy(particleSource->ctx);
    free(particleSource);
    blParticleSourceDestroy(next);
  }
}

/*
 * Implementation of the uniform particle source
 */
struct UniformCtx {
  struct BBox volume;
  int numPtcls;
  double vbar[3];
  double deltaV[3];
  int internalStateSize;
  double *internalState;
};

static int uniformGetNumParticles(void *ctx) {
  struct UniformCtx *uctx = ctx;
  return uctx->numPtcls;
}

static void uniformCreateParticles(double *x, double *y, double *z,
      double *vx, double *vy, double *vz, double *internalState, void *ctx) {
  struct UniformCtx *uCtx = ctx;
  const struct BBox *box = &uCtx->volume;
  int i;

  for (i = 0; i < uCtx->numPtcls; ++i) {
    x[i] = box->xmin + (box->xmax - box->xmin) * sprng();
    y[i] = box->ymin + (box->ymax - box->ymin) * sprng();
    z[i] = box->zmin + (box->zmax - box->zmin) * sprng();
    vx[i] = generateGaussianNoise(uCtx->vbar[0], uCtx->deltaV[0]);
    vy[i] = generateGaussianNoise(uCtx->vbar[1], uCtx->deltaV[1]);
    vz[i] = generateGaussianNoise(uCtx->vbar[2], uCtx->deltaV[2]);
    memcpy(&internalState[i * uCtx->internalStateSize],
        uCtx->internalState, uCtx->internalStateSize * sizeof(double));
  }
}

void uniformDestroy(void *ctx) {
  struct UniformCtx *uCtx = ctx;
  free(uCtx->internalState);
  free(uCtx);
}

struct ParticleSource *blParticleSourceUniformCreate(
    struct BBox volume, int numPtcls, double *vbar, double *deltaV,
    int internalStateSize, double *internalState,
    struct ParticleSource *next) {
  struct ParticleSource *self;
  blParticleSourceCreate(next, &self);
  self->getNumParticles = uniformGetNumParticles;
  self->createParticles = uniformCreateParticles;
  self->destroy = uniformDestroy;
  struct UniformCtx *ctx = malloc(sizeof(struct UniformCtx));
  ctx->volume = volume;
  ctx->numPtcls = numPtcls;
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

