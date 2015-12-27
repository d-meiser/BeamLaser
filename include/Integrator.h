#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Errors.h>


typedef struct BLIntegrator_ *BLIntegrator;

typedef void (*BLIntegratorRHS)(double t, const double *x, double *y,
    void *ctx);

BL_STATUS blIntegratorCreate(const char* name, int n, BLIntegrator *integrator);
void blIntegratorDestroy(BLIntegrator *integrator);
void blIntegratorTakeStep(BLIntegrator integrator, double t, double dt,
    BLIntegratorRHS rhs, void *ctx);

#endif
