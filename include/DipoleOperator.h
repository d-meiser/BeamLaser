#ifndef DIPOLE_OPERATOR_H
#define DIPOLE_OPERATOR_H

struct DipoleOperator {
  void (*apply)(int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double psi, double *result, double *polarization);
  void *ctx;
};

void blDipoleOperatorApply(struct DipoleOperator *op,
    int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double psi, double *result, double *polarization);

void blDipoleOperatorDestroy(struct DipoleOperator *op);

#endif

