#include <Integrator.h>

BL_STATUS blIntegratorCreate(const char* name, int n,
    BLIntegrator *integrator) {
  return BL_UNKNOWN_INTEGRATOR;
}

void blIntegratorDestroy(BLIntegrator *integrator);
void blIntegratorTakeStep(BLIntegrator integrator, double t, double dt,
    BLIntegratorRHS rhs, void *ctx);
