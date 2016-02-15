#include <FieldUpdate.h>
#include <stdlib.h>
#include <math.h>


struct FieldUpdateCtx {
  double damping;
  double noise;
};

static void fieldUpdateTakeStep(
    double t, double dt, struct BLSimulationState *state, void *c) {
  BL_UNUSED(t);
  struct FieldUpdateCtx *ctx = (struct FieldUpdateCtx*)c;
  
  state->fieldState.q *= exp(-0.5 * ctx->damping * dt);
  state->fieldState.p *= exp(-0.5 * ctx->damping * dt);

  state->fieldState.q += blGenerateGaussianNoise(0.0, sqrt(dt * ctx->noise));
  state->fieldState.p += blGenerateGaussianNoise(0.0, sqrt(dt * ctx->noise));

  state->fieldState.q *= exp(-0.5 * ctx->damping * dt);
  state->fieldState.p *= exp(-0.5 * ctx->damping * dt);
}

static void fieldUpdateDestroy(void *c) {
  free(c);
}

struct BLUpdate *blFieldUpdateCreate(double damping, double noise) {
  struct BLUpdate *this = malloc(sizeof(*this));
  this->takeStep = fieldUpdateTakeStep;
  this->destroy = fieldUpdateDestroy;
  struct FieldUpdateCtx *ctx = malloc(sizeof(*ctx));
  ctx->damping = damping;
  ctx->noise = noise;
  this->ctx = ctx;
  return this;
}

