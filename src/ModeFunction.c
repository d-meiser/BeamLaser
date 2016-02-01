#include <ModeFunction.h>
#include <stdlib.h>
#include <math.h>
#include <Utilities.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define H_BAR 1.0545718e-34
#define EPSILON_0 8.85e-12
#define SPEED_OF_LIGHT 299792458.0

#define MODE_FUN_CHUNK_SIZE 16


void blModeFunctionDestroy(struct BLModeFunction *modeFunction) {
  modeFunction->destroy(modeFunction->ctx);
  free(modeFunction);
}

void blModeFunctionEvaluate(struct BLModeFunction *modeFunction,
    int n,
    const double * restrict x,
    const double * restrict y,
    const double * restrict z,
    double complex * restrict fx,
    double complex * restrict fy,
    double complex * restrict fz) {
  modeFunction->evaluate(n, x, y, z, fx, fy, fz, modeFunction->ctx);
}

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
    void *c) {
  struct ModeFunctionSimplifiedGaussianCtx *ctx =
    (struct ModeFunctionSimplifiedGaussianCtx*)c;
  double waist = ctx->waist;
  double k = ctx->k;
  double amplitude = ctx->amplitude;

  int i, j;
  int numChunks = n / MODE_FUN_CHUNK_SIZE;
  for (i = 0; i < numChunks; ++i) {

    const double * restrict xx = x + i * MODE_FUN_CHUNK_SIZE;
    const double * restrict yy = y + i * MODE_FUN_CHUNK_SIZE;
    const double * restrict zz = z + i * MODE_FUN_CHUNK_SIZE;

    double gaussianTerm[MODE_FUN_CHUNK_SIZE];
    for (j = 0; j < MODE_FUN_CHUNK_SIZE; ++j) {
      gaussianTerm[j] = -0.5 * (yy[j] * yy[j] + zz[j] * zz[j]) / (waist * waist);
    }
    for (j = 0; j < MODE_FUN_CHUNK_SIZE; ++j) {
      gaussianTerm[j] = exp(gaussianTerm[j]);
    }

    double sineTerm[MODE_FUN_CHUNK_SIZE];
    for (j = 0; j < MODE_FUN_CHUNK_SIZE; ++j) {
      sineTerm[j] = k * xx[j];
    }
    for (j = 0; j < MODE_FUN_CHUNK_SIZE; ++j) {
      sineTerm[j] = sin(sineTerm[j]);
    }

    double complex * restrict fyy = fy + i * MODE_FUN_CHUNK_SIZE;
    for (j = 0; j < MODE_FUN_CHUNK_SIZE; ++j) {
      fyy[j] = amplitude * gaussianTerm[j] * sineTerm[j];
    }
  }

  /* clean up of peeled loop */
  for (i = numChunks * MODE_FUN_CHUNK_SIZE; i < n; ++i) {
    double arge = -0.5 *(y[i] * y[i] + z[i] * z[i]) / (waist * waist);
    double args = k * x[i];
    fy[i] = amplitude * exp(arge) * sin(args);
  }

  for (i = 0; i < n; ++i) {
    fx[i] = 0;
  }

  for (i = 0; i < n; ++i) {
    fz[i] = 0;
  }
}

static void modeFunctionSimplifiedGaussianDestroy(void *ctx) {
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
  ctx->amplitude = sqrt(omega / (2.0 * H_BAR * EPSILON_0 * veff));
  return this;
}


struct ModeFunctionUniformCtx {
  double fx;
  double fy;
  double fz;
};

static void modeFunctionUniformEvaluate(int n,
    const double * restrict x,
    const double * restrict y,
    const double * restrict z,
    double complex * restrict fx,
    double complex * restrict fy,
    double complex * restrict fz,
    void *c) {
  BL_UNUSED(x);
  BL_UNUSED(y);
  BL_UNUSED(z);
  struct ModeFunctionUniformCtx *ctx = (struct ModeFunctionUniformCtx*)c;
  int i;
  for (i = 0; i < n; ++i) {
    fx[i] = ctx->fx;
  }
  for (i = 0; i < n; ++i) {
    fy[i] = ctx->fy;
  }
  for (i = 0; i < n; ++i) {
    fz[i] = ctx->fz;
  }
}

static void modeFunctionUniformDestroy(void *ctx) {
  free(ctx);
}

struct BLModeFunction *blModeFunctionUniformCreate(double fx, double fy, double fz) {
  struct BLModeFunction *this = malloc(sizeof(*this));
  this->evaluate = modeFunctionUniformEvaluate;
  this->destroy = modeFunctionUniformDestroy;
  struct ModeFunctionUniformCtx *ctx = malloc(sizeof(*ctx));
  this->ctx = ctx;
  ctx->fx = fx;
  ctx->fy = fy;
  ctx->fz = fz;
  return this;
}
