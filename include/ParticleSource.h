#ifndef PARTICLE_SOURCE_H
#define PARTICLE_SOURCE_H

#include <Utilities.h>


struct ParticleSource {
  int (*getNumParticles)(void *ctx);
  void (*createParticles)(double *x, double *y, double *z,
      double *vx, double *vy, double *vz, double *internalState, void *ctx);
  void (*destroy)(void *ctx);
  void *ctx;
  struct ParticleSource *next;
};

int blParticleSourceGetNumParticles(struct ParticleSource *particleSource);
void blParticleSourceCreateParticles(struct ParticleSource *particleSource,
    double *x, double *y, double *z, double *vx, double *vy, double *vz,
    double *internalState);
void blParticleSourceDestroy(struct ParticleSource *particleSource);

struct ParticleSource *blParticleSourceUniformCreate(
    struct BBox volume, int nbar,
    double *vbar, double *deltaV,
    int internalStateSize, double *internalState,
    struct ParticleSource *next);

#endif

