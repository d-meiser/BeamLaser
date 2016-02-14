/*
Copyright 2014 Dominic Meiser

This file is part of BeamLaser.

BeamLaser is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

BeamLaser is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with BeamLaser.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <Integrator.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static BL_STATUS integratorCreateRK4(BLIntegrator integrator);
static void integratorDestroyRK4(BLIntegrator integrator);
static void integratorTakeStepRK4(BLIntegrator integrator, double t, double dt, int n,
    BLIntegratorRHS rhs, const double *x, double *y, void *ctx);

struct BLIntegrator_ {
  int n;
  void (*destroy)(struct BLIntegrator_ *this);
  void (*takeStep)(struct BLIntegrator_ *this, double t, double dt, int n,
        BLIntegratorRHS rhs, const double *x, double *y, void *ctx);
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

void blIntegratorTakeStep(BLIntegrator integrator, double t, double dt, int n,
    BLIntegratorRHS rhs, const double *x, double *y, void *ctx) {
  integrator->takeStep(integrator, t, dt, n, rhs, x, y, ctx);
}


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

static void daxpy(double * restrict w, double alpha,
    const double * restrict x, const double * restrict y,
    int dim)
{
  int i;
  for (i = 0; i < dim; ++i) {
    w[i] = alpha * x[i] + y[i];
  }
}

static void integratorTakeStepRK4(BLIntegrator integrator, double t, double dt,
                                  int n, BLIntegratorRHS rhs,
                                  const double *x, double *y, void *ctx) {
  struct RK4Ctx *rk4Ctx = integrator->intCtx;
  double prefactor;
  int i;
  double * restrict k1 = rk4Ctx->k1;
  double * restrict k2 = rk4Ctx->k2;
  double * restrict k3 = rk4Ctx->k3;
  double * restrict k4 = rk4Ctx->k4;

  rhs(t, n, x, k1, ctx);
  daxpy(rk4Ctx->tmp, 0.5 * dt, k1, x, n);
  rhs(t + 0.5 * dt, n, rk4Ctx->tmp, k2, ctx);
  daxpy(rk4Ctx->tmp, 0.5 * dt, k2, x, n);
  rhs(t + 0.5 * dt, n, rk4Ctx->tmp, k3, ctx);
  daxpy(rk4Ctx->tmp, dt, k3, x, n);
  rhs(t + dt, n, rk4Ctx->tmp, k4, ctx);
  prefactor = dt / 6.0;
  if (y == x) {
    double * restrict yp = y;
    for (i = 0; i < n; ++i) {
      yp[i] += prefactor * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]);
    }
  } else {
    const double * restrict xp = x;
    double * restrict yp = y;
    for (i = 0; i < n; ++i) {
      yp[i] = xp[i] + prefactor * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]);
    }
  }
}
