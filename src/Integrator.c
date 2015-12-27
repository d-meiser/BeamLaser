#include <Integrator.h>
#include <stdlib.h>
#include <string.h>

static BL_STATUS integratorCreateRK4(BLIntegrator integrator);
static void integratorDestroyRK4(BLIntegrator integrator);
static void integratorTakeStepRK4(BLIntegrator integrator, double t, double dt,
    BLIntegratorRHS rhs, void *ctx);

struct BLIntegrator_ {
  int n;
  void (*destroy)(struct BLIntegrator_ *this);
  void (*takeStep)(struct BLIntegrator_ *this, double t, double dt,
        BLIntegratorRHS rhs, void *ctx);
  void *intCtx;
};

BL_STATUS blIntegratorCreate(const char* name, int n,
    BLIntegrator *integrator) {
  *integrator = malloc(sizeof(**integrator));
  (*integrator)->n = n;
  (*integrator)->intCtx = 0;
  if (!strcmp(name, "RK4")) {
    return integratorCreateRK4(*integrator);
  } else {
    return BL_UNKNOWN_INTEGRATOR;
  }
}

void blIntegratorDestroy(BLIntegrator *integrator) {
  (*integrator)->destroy(*integrator);
  free(*integrator);
  *integrator = 0;
}

void blIntegratorTakeStep(BLIntegrator integrator, double t, double dt,
    BLIntegratorRHS rhs, void *ctx);


/*
 * Implementation for RK4
 **/
struct RK4Ctx {
  double *k1, *k2, *k3, *k4, *tmp;
};

BL_STATUS integratorCreateRK4(struct BLIntegrator_ *integrator) {
  struct RK4Ctx *ctx = malloc(sizeof(*ctx));
  if (!ctx) return BL_OUT_OF_MEMORY;
  ctx->k1 = malloc(integrator->n * sizeof(*ctx->k1));
  if (!ctx->k1) return BL_OUT_OF_MEMORY;
  ctx->k2 = malloc(integrator->n * sizeof(*ctx->k2));
  if (!ctx->k2) return BL_OUT_OF_MEMORY;
  ctx->k3 = malloc(integrator->n * sizeof(*ctx->k3));
  if (!ctx->k3) return BL_OUT_OF_MEMORY;
  ctx->k4 = malloc(integrator->n * sizeof(*ctx->k4));
  if (!ctx->k4) return BL_OUT_OF_MEMORY;
  ctx->tmp = malloc(integrator->n * sizeof(*ctx->tmp));
  if (!ctx->tmp) return BL_OUT_OF_MEMORY;
  integrator->intCtx = ctx;

  integrator->destroy = integratorDestroyRK4;
  integrator->takeStep = integratorTakeStepRK4;

  return BL_SUCCESS;
}

static void integratorDestroyRK4(BLIntegrator integrator) {
  struct RK4Ctx *ctx = integrator->intCtx;
  free(ctx->k1);
  free(ctx->k2);
  free(ctx->k3);
  free(ctx->k4);
  free(ctx->tmp);
  free(ctx);
  integrator->intCtx = 0;
}

static void integratorTakeStepRK4(BLIntegrator integrator, double t, double dt,
    BLIntegratorRHS rhs, void *ctx) {
}
