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

void blDipoleOperatorApply(struct BLDipoleOperator *op,
    int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result, double *polarization) {
  op->apply(stride, numPtcls, ex, ey, ez, psi, result, polarization, op->ctx);
}

void blDipoleOperatorApplyNoPolarization(struct BLDipoleOperator *op,
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

void blDipoleOperatorDestroy(struct BLDipoleOperator *op) {
  op->destroy(op->ctx);
  free(op);
}

static void dipoleOperatorTLAApply(int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result, double *polarization,
    void *ctx) {
  BL_UNUSED(ctx);
  BL_UNUSED(ez);
  double complex y[2];
  double complex pol = 0;
  double complex *x, *r;
  int i, j;
  for (i = 0; i < numPtcls; ++i) {
    double complex ep = (ex[i] + I * ey[i]) * M_SQRT1_2;
    x = (double complex*)(psi + i * stride);
    r = (double complex*)(result + i * stride);
    for (j = 0; j < 2; ++j) {
      y[j] = 0;
    }
    y[0] = ep * x[1];
    y[1] = conj(ep) * x[0];
    for (j = 0; j < 2; ++j) {
      pol += conj(x[j]) * y[j];
    }
    for (j = 0; j < 2; ++j) {
      r[j] = y[j];
    }
  }
  *((double complex*)polarization) = pol;
}

static void dipoleOperatorTLAApplyNoPolarization(int stride, int numPtcls,
                                   const double *ex, const double *ey, const double *ez,
                                   const double *psi, double *result,
                                   void *ctx) {
  BL_UNUSED(ctx);
  BL_UNUSED(ez);
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
    y[0] = ep * x[1];
    y[1] = conj(ep) * x[0];
    for (j = 0; j < 2; ++j) {
      r[j] = y[j];
    }
  }
}

static void dipoleOperatorTLADestroy(void *ctx) {
  BL_UNUSED(ctx);
}

struct BLDipoleOperator *blDipoleOperatorTLACreate() {
  struct BLDipoleOperator *op = malloc(sizeof(*op));
  op->apply = dipoleOperatorTLAApply;
  op->applyNoPolarization = dipoleOperatorTLAApplyNoPolarization;
  op->destroy = dipoleOperatorTLADestroy;
  op->ctx = 0;
  return op;
}
