#include <FieldUpdate.h>
#include <stdlib.h>
#include <math.h>


struct FieldUpdateCtx {
  double detuning;
  double damping;
  double noise;
};

static void fieldUpdateTakeStep(
    double t, double dt, struct BLSimulationState *state, void *ct) {
  BL_UNUSED(t);
  struct FieldUpdateCtx *ctx = (struct FieldUpdateCtx*)ct;
  
  double c = cos(0.5 * ctx->detuning * dt);
  double s = sin(0.5 * ctx->detuning * dt);
  double e = exp(-0.5 * ctx->damping * dt);

  double qTmp;
  double pTmp;
  qTmp = e * (c * state->fieldState.q - s * state->fieldState.p);
  pTmp = e * (c * state->fieldState.p + s * state->fieldState.q);

  qTmp += blGenerateGaussianNoise(0.0, sqrt(dt * ctx->noise));
  pTmp += blGenerateGaussianNoise(0.0, sqrt(dt * ctx->noise));

  state->fieldState.q = e * (c * qTmp - s * pTmp);
  state->fieldState.p = e * (c * pTmp + s * qTmp);
}

static void fieldUpdateDestroy(void *c) {
  free(c);
}

struct BLUpdate *blFieldUpdateCreate(double detuning, double damping,
    double noise) {
  struct BLUpdate *this = malloc(sizeof(*this));
  this->takeStep = fieldUpdateTakeStep;
  this->destroy = fieldUpdateDestroy;
  struct FieldUpdateCtx *ctx = malloc(sizeof(*ctx));
  ctx->detuning = detuning;
  ctx->damping = damping;
  ctx->noise = noise;
  this->ctx = ctx;
  return this;
}

