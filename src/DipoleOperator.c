#include <DipoleOperator.h>
#include <Utilities.h>
#include <stdlib.h>

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
