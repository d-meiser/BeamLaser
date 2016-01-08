#include <DipoleOperator.h>
#include <Utilities.h>
#include <stdlib.h>
#include <complex.h>

#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif

void blDipoleOperatorApply(struct DipoleOperator *op,
    int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result, double *polarization) {
  op->apply(stride, numPtcls, ex, ey, ez, psi, result, polarization, op->ctx);
}

void blDipoleOperatorApplyNoPolarization(struct DipoleOperator *op,
    int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result) {
  if (op->applyNoPolarization) {
    op->applyNoPolarization(stride, numPtcls, ex, ey, ez, psi, result, op->ctx);
  } else { 
    double polarization[2];
    op->apply(stride, numPtcls, ex, ey, ez, psi, result, polarization, op->ctx);
  }
}

void blDipoleOperatorDestroy(struct DipoleOperator *op) {
  op->destroy(op->ctx);
  free(op);
}

static void dipoleOperatorTLAApply(int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result, double *polarization,
    void *ctx) {
  double complex y[2];
  double complex *x, *r;
  int i, j;
  for (i = 0; i < numPtcls; ++i) {
    double complex ep = (ex[i] + I * ey[i]) * M_SQRT1_2;
    x = (double complex*)(psi + i * stride);
    r = (double complex*)(result + i * stride);
    for (j = 0; j < 2; ++j) {
      y[j] = 0;
    }
    for (j = 0; j < 2; ++j) {
      r[j] = y[j];
    }
  }
  BL_UNUSED(stride);
  BL_UNUSED(numPtcls);
  BL_UNUSED(ex);
  BL_UNUSED(ey);
  BL_UNUSED(ez);
  BL_UNUSED(psi);
  BL_UNUSED(result);
  BL_UNUSED(polarization);
  BL_UNUSED(ctx);
}

static void dipoleOperatorTLADestroy(void *ctx) {
  BL_UNUSED(ctx);
}

struct DipoleOperator *blDipoleOperatorTLACreate() {
  struct DipoleOperator *op = malloc(sizeof(*op));
  op->apply = dipoleOperatorTLAApply;
  op->applyNoPolarization = 0;
  op->destroy = dipoleOperatorTLADestroy;
  op->ctx = 0;
  return op;
}
