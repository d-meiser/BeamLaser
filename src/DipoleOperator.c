#include <DipoleOperator.h>
#include <Utilities.h>
#include <stdlib.h>
#include <complex.h>


void blDipoleOperatorApply(struct BLDipoleOperator *op,
    int stride, int numPtcls,
    const double complex *ex,
    const double complex *ey,
    const double complex *ez,
    const double complex *psi, double complex *result) {
  op->apply(stride, numPtcls, ex, ey, ez, psi, result, op->ctx);
}

void blDipoleOperatorComputeD(struct BLDipoleOperator *op,
    int stride, int numPtcls,
    const double complex *psi,
    double complex *dx, double complex *dy, double complex *dz) {
  op->computeD(stride, numPtcls, psi, dx, dy, dz, op->ctx);
}

void blDipoleOperatorDestroy(struct BLDipoleOperator *op) {
  op->destroy(op->ctx);
  free(op);
}


struct TLACtx {
  double dipoleMatrixElement;
};

static void dipoleOperatorTLAApply(int stride, int numPtcls,
    const double complex *ex,
    const double complex *ey,
    const double complex *ez,
    const double complex *psi, double complex *result,
    void *c) {
  BL_UNUSED(ex);
  BL_UNUSED(ez);
  struct TLACtx *ctx = c;
  double dipoleMatrixElement = ctx->dipoleMatrixElement;
  double complex y[2];
  const double complex *x;
  double complex *r;
  int i;
  for (i = 0; i < numPtcls; ++i) {
    x = psi + i * stride;
    r = result + i * stride;

    y[0] = 0;
    y[1] = 0;

    double complex ep = dipoleMatrixElement * ey[i];
    y[0] = conj(ep) * x[1];
    y[1] = ep * x[0];

    r[0] = y[0];
    r[1] = y[1];
  }
}

static void dipoleOperatorTLAComputeD(int stride, int numPtcls,
    const double complex *psi,
    double complex *dx, double complex *dy, double complex *dz, void *c) {
  BL_UNUSED(dz);
  BL_UNUSED(dy);
  struct TLACtx *ctx = c;
  double dipoleMatrixElement = ctx->dipoleMatrixElement;
  const double complex *x;
  int i, j;
  for (i = 0; i < numPtcls; ++i) {
    dx[i] = 0;
  }
  for (i = 0; i < numPtcls; ++i) {
    x = psi + i * stride;
    const double *xr = (const double*)x;
    double nrmSquared = 0;
    for (j = 0; j < 4; ++j) {
      nrmSquared += xr[j] * xr[j];
    }
    dy[i] = conj(x[0]) * dipoleMatrixElement * x[1] / nrmSquared;
  }
  for (i = 0; i < numPtcls; ++i) {
    dz[i] = 0;
  }
}

static void dipoleOperatorTLADestroy(void *ctx) {
  free(ctx);
}

struct BLDipoleOperator *blDipoleOperatorTLACreate(double dipoleMatrixElement) {
  struct BLDipoleOperator *op = malloc(sizeof(*op));
  op->apply = dipoleOperatorTLAApply;
  op->computeD = dipoleOperatorTLAComputeD;
  op->destroy = dipoleOperatorTLADestroy;
  struct TLACtx *ctx = malloc(sizeof(*ctx));
  ctx->dipoleMatrixElement = dipoleMatrixElement;
  op->ctx = ctx;
  return op;
}
