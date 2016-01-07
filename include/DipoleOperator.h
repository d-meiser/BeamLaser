#ifndef DIPOLE_OPERATOR_H
#define DIPOLE_OPERATOR_H

struct DipoleOperator {
  void (*apply)(int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result, double *polarization,
    void *ctx);
  void (*applyNoPolarization)(int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result,
    void *ctx);
  void (*destroy)(void *ctx);
  void *ctx;
};

void blDipoleOperatorApply(struct DipoleOperator *op,
    int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result, double *polarization);

void blDipoleOperatorApplyNoPolarization(struct DipoleOperator *op,
    int stride, int numPtcls,
    const double *ex, const double *ey, const double *ez,
    const double *psi, double *result);

void blDipoleOperatorDestroy(struct DipoleOperator *op);


struct DipoleOperator *blDipoleOperatorTLACreate();

#endif

