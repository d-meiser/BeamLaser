#include <ModeFunction.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define H_BAR 1.0545718e-34
#define EPSILON_0 8.85e-12
#define SPEED_OF_LIGHT 299792458.0


struct ModeFunctionSimplifiedGaussianCtx {
  double waist;
  double k;
  double amplitude;
};

static void modeFunctionSimplifiedGaussianEvaluate(int n,
    const double * restrict x,
    const double * restrict y,
    const double * restrict z,
    double complex * restrict fx,
    double complex * restrict fy,
    double complex * restrict fz,
    void *ctx) {
}

void modeFunctionSimplifiedGaussianDestroy(void *ctx) {
  free(ctx);
}

struct BLModeFunction *blModeFunctionSimplifiedGaussianCreate(
    double waist, double lambda, double length) {
  struct BLModeFunction *this = malloc(sizeof(*this));
  this->evaluate = modeFunctionSimplifiedGaussianEvaluate;
  this->destroy = modeFunctionSimplifiedGaussianDestroy;
  struct ModeFunctionSimplifiedGaussianCtx *ctx = malloc(sizeof(*ctx));
  this->ctx = ctx;
  ctx->waist = waist;
  ctx->k = 2.0 * M_PI / lambda;
  double veff = waist * waist * length;
  double omega = ctx->k * SPEED_OF_LIGHT;
  ctx->amplitude = sqrt(H_BAR * omega / (2.0 * EPSILON_0 * veff));
  return this;
}
