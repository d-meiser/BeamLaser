#ifndef MODE_FUNCTION_H
#define MODE_FUNCTION_H

#include <complex.h>

struct BLModeFunction {
  void (*evaluate)(int n,
    const double * restrict x,
    const double * restrict y,
    const double * restrict z,
    double complex * restrict fx,
    double complex * restrict fy,
    double complex * restrict fz,
    void *ctx);
  void (*destroy)(void *);
  void *ctx;
};

void blModeFunctionEvaluate(struct BLModeFunction *modeFunction,
    int n,
    const double * restrict x,
    const double * restrict y,
    const double * restrict z,
    double complex * restrict fx,
    double complex * restrict fy,
    double complex * restrict fz);
void blModeFunctionDestroy(struct BLModeFunction *modeFunction);

struct BLModeFunction *blModeFunctionSimplifiedGaussianCreate(
    double waist, double lambda, double length);

struct BLModeFunction *blModeFunctionUniformCreate(double fx, double fy, double fz);

#endif

