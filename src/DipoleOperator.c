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

void blDipoleOperatorExpectationValue(struct BLDipoleOperator *op,
    int stride, int numPtcls,
    const double complex *psi,
    double complex *dx, double complex *dy, double complex *dz) {
  op->expectationValue(stride, numPtcls, psi, dx, dy, dz, op->ctx);
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
  BL_UNUSED(ez);
  BL_UNUSED(ey);
  struct TLACtx *ctx = c;
  double dipoleMatrixElement = ctx->dipoleMatrixElement;
  double complex y[2];
  const double complex *x;
  double complex *r;
  int i, j;
  for (i = 0; i < numPtcls; ++i) {
    x = psi + i * stride;
    r = result + i * stride;

    for (j = 0; j < 2; ++j) {
      y[j] = 0;
    }

    double complex ep = dipoleMatrixElement * ex[i];
    y[0] = conj(ep) * x[1];
    y[1] = ep * x[0];

    for (j = 0; j < 2; ++j) {
      r[j] = y[j];
    }
  }
}

static void dipoleOperatorTLADestroy(void *ctx) {
  free(ctx);
}

struct BLDipoleOperator *blDipoleOperatorTLACreate(double dipoleMatrixElement) {
  struct BLDipoleOperator *op = malloc(sizeof(*op));
  op->apply = dipoleOperatorTLAApply;
  op->destroy = dipoleOperatorTLADestroy;
  struct TLACtx *ctx = malloc(sizeof(*ctx));
  ctx->dipoleMatrixElement = dipoleMatrixElement;
  op->ctx = ctx;
  return op;
}
