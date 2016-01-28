#ifndef DIPOLE_OPERATOR_H
#define DIPOLE_OPERATOR_H

#include <complex.h>


struct BLDipoleOperator {
  void (*apply)(int stride, int numPtcls,
      const double complex *ex,
      const double complex *ey,
      const double complex *ez,
      const double complex *psi, double complex *result, void *ctx);
  void (*computeD)(int stride, int numPtcls,
      const double complex *psi,
      double complex *dx, double complex *dy, double complex *dz,
      void *ctx);
  void (*destroy)(void *ctx);
  void *ctx;
};

void blDipoleOperatorApply(struct BLDipoleOperator *op,
    int stride, int numPtcls,
    const double complex *ex,
    const double complex *ey,
    const double complex *ez,
    const double complex *psi, double complex *result);

void blDipoleOperatorComputeD(struct BLDipoleOperator *op,
    int stride, int numPtcls, const double complex *psi,
    double complex *dx, double complex *dy, double complex *dz);

void blDipoleOperatorExpectationValue(struct BLDipoleOperator *op,
    int stride, int numPtcls,
    const double complex *psi,
    double complex *dx,
    double complex *dy,
    double complex *dz);

void blDipoleOperatorDestroy(struct BLDipoleOperator *op);


struct BLDipoleOperator *blDipoleOperatorTLACreate(double dipoleMatrixElement);

#endif

