#include <BeamLaser.h>
#include <config.h>

#ifdef BL_WITH_MPI
#include <mpi.h>
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifndef BL_WITH_MPI
#define MPI_Request int
#endif

#define BL_UNUSED(a) (void)(a)

static const int INTERNAL_STATE_DIM = 4;

struct FieldState {
  double q;
  double p;
};

struct SimulationState {
  double t;
  struct FieldState fieldState;
  struct BLEnsemble ensemble;
};

struct Configuration {
  int numSteps;
  double particleWeight;
  double dipoleMatrixElement;
  double nbar;
  int maxNumParticles;
  double dt;
  double vbar;
  double deltaV;
  double alpha;
  double kappa;
  struct BBox simulationDomain;
  struct BBox sourceVolume;
};

struct IntegratorCtx {
  struct BLEnsemble *ensemble;
  double *x;
  double *ex;
  double *ey;
  double *ez;
};

void setDefaults(struct Configuration *conf);
void computeSourceVolume(struct Configuration *conf);
void particleSink(const struct Configuration *conf, struct BLEnsemble *ensemble);
void particleSource(const struct Configuration *conf, struct BLEnsemble *ensemble);
void blFieldUpdate(double dt, double kappa, struct FieldState *fieldState);
void blFieldAtomInteraction(double dt, struct FieldState *fieldState,
        struct IntegratorCtx *integratorCtx, BLIntegrator integrator);
void scatterFieldEnd(MPI_Request req, const struct FieldState *fieldState,
    double *fieldDest);
MPI_Request scatterFieldBegin(const struct FieldState *fieldState,
    double *fieldDest);
void interactionRHS(double t, int n, const double *x, double *y,
                    void *ctx);
static void modeFunction(double x, double y, double z,
                         double *fx, double *fy, double *fz);

int main(int argn, char **argv) {
#ifdef BL_WITH_MPI
  MPI_Init(&argn, &argv);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif
  struct SimulationState simulationState;
  struct Configuration conf;
  struct IntegratorCtx integratorCtx;
  BLIntegrator integrator;
  BL_STATUS stat;
  int i;

  setDefaults(&conf);
  computeSourceVolume(&conf);
  blIntegratorCreate("RK4", conf.maxNumParticles * INTERNAL_STATE_DIM,
                     &integrator);
  
  stat = blEnsembleInitialize(conf.maxNumParticles, INTERNAL_STATE_DIM,
      &simulationState.ensemble);
  integratorCtx.ensemble = &simulationState.ensemble;
  integratorCtx.x = malloc(
      INTERNAL_STATE_DIM * conf.maxNumParticles * sizeof(double));
  integratorCtx.ex = malloc(conf.maxNumParticles * sizeof(double));
  integratorCtx.ey = malloc(conf.maxNumParticles * sizeof(double));
  integratorCtx.ez = malloc(conf.maxNumParticles * sizeof(double));

  if (stat != BL_SUCCESS) return stat;
  simulationState.fieldState.q = 1.0;
  simulationState.fieldState.p = 0.0;
  if (stat != BL_SUCCESS) return stat;

  for (i = 0; i < conf.numSteps; ++i) {
    particleSink(&conf, &simulationState.ensemble);
    particleSource(&conf, &simulationState.ensemble);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blFieldAtomInteraction(conf.dt, &simulationState.fieldState,
        &integratorCtx, integrator);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    printf("%d %le %le\n", i, simulationState.fieldState.q, simulationState.fieldState.p);
  }

  free(integratorCtx.x);
  free(integratorCtx.ex);
  free(integratorCtx.ey);
  free(integratorCtx.ez);
  blIntegratorDestroy(&integrator);
  blEnsembleFree(&simulationState.ensemble);

#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return BL_SUCCESS;
}

void setDefaults(struct Configuration *conf) {
  conf->numSteps = 10;
  conf->particleWeight = 1.0e6;
  conf->dipoleMatrixElement = 1.0e-5 * 1.0e-29;
  conf->nbar = 1.0e3;
  conf->maxNumParticles = 2000;
  conf->dt = 1.0e-8;
  conf->vbar = 3.0e2;
  conf->deltaV = 1.0e1;
  conf->alpha = 1.0e-2;
  conf->kappa = 2.0 * M_PI * 1.0e6;
  conf->simulationDomain.xmin = -1.0e-4;
  conf->simulationDomain.xmax = 1.0e-4;
  conf->simulationDomain.ymin = -1.0e-4;
  conf->simulationDomain.ymax = 1.0e-4;
  conf->simulationDomain.zmin = -1.0e-4;
  conf->simulationDomain.zmax = 1.0e-4;
}

void computeSourceVolume(struct Configuration *conf) {
  conf->sourceVolume.xmin = conf->simulationDomain.xmin;
  conf->sourceVolume.xmax = conf->simulationDomain.xmax;
  conf->sourceVolume.ymin = conf->simulationDomain.ymin;
  conf->sourceVolume.ymax = conf->simulationDomain.ymax;
  conf->sourceVolume.zmin = conf->simulationDomain.zmax;
  conf->sourceVolume.zmax = conf->simulationDomain.zmax + conf->vbar * conf->dt;
}

void particleSink(const struct Configuration *conf,
                  struct BLEnsemble *ensemble) {
  blEnsembleRemoveBelow(conf->simulationDomain.zmin, ensemble->z, ensemble);
}

void particleSource(const struct Configuration *conf,
                    struct BLEnsemble *ensemble) {
  int numCreate = round(
      conf->nbar * 
      (conf->dt * conf->vbar) /
      (conf->simulationDomain.zmax - conf->simulationDomain.zmin)
      );
  int i, j;
  while (numCreate > 0) {
    i = blRingBufferAppendOne(&ensemble->buffer);
    blEnsembleCreateParticle(conf->sourceVolume, conf->vbar, conf->deltaV,
                             conf->alpha, i, ensemble);
    for (j = 0; j < ensemble->internalStateSize; ++j) {
      ensemble->internalState[i * ensemble->internalStateSize + j] = 0;
    }
    ensemble->internalState[i * ensemble->internalStateSize + 2] = 1.0;
    --numCreate;
  }
}

void blFieldUpdate(double dt, double kappa, struct FieldState *fieldState) {
  fieldState->q = fieldState->q * exp(-0.5 * kappa * dt);
  fieldState->p = fieldState->p * exp(-0.5 * kappa * dt);
}

void blFieldAtomInteraction(double dt, struct FieldState *fieldState,
        struct IntegratorCtx *integratorCtx, BLIntegrator integrator) {
  int i, ip;
  struct BLEnsemble *ensemble = integratorCtx->ensemble;
  int n = blRingBufferSize(ensemble->buffer) * ensemble->internalStateSize + 2;
  double *x = integratorCtx->x;

  /* 
   * Pack field and internal state data in contiguous buffer;
   * Field is replicated and integrated redundantly
   */
  MPI_Request fieldRequest = scatterFieldBegin(fieldState, x);
  for (i = 0, ip = ensemble->buffer.begin; ip != ensemble->buffer.end;
      ++i, ip = blRingBufferNext(ensemble->buffer, ip)) {
    memcpy(&x[2 + i * INTERNAL_STATE_DIM],
        &ensemble->internalState[ip * INTERNAL_STATE_DIM],
        INTERNAL_STATE_DIM * sizeof(double));
  }
  for (i = 0, ip = ensemble->buffer.begin; ip != ensemble->buffer.end;
      ++i, ip = blRingBufferNext(ensemble->buffer, ip)) {
    modeFunction(ensemble->x[ip], ensemble->y[ip], ensemble->z[ip],
        &integratorCtx->ex[i], &integratorCtx->ey[i], &integratorCtx->ez[i]);
  }

  scatterFieldEnd(fieldRequest, fieldState, x);

  blIntegratorTakeStep(integrator, 0.0, dt, n, interactionRHS, x, x,
      integratorCtx);

  fieldState->q = x[0];
  fieldState->p = x[1];
  for (i = 0, ip = ensemble->buffer.begin; ip != ensemble->buffer.end;
      ++i, ip = blRingBufferNext(ensemble->buffer, ip)) {
    memcpy(&ensemble->internalState[ip * ensemble->internalStateSize],
           &x[2 + i * ensemble->internalStateSize],
           ensemble->internalStateSize * sizeof(double));
  }
}

void interactionRHS(double t, int n, const double *x, double *y,
                    void *ctx) {
  BL_UNUSED(t);
  BL_UNUSED(n);
  struct IntegratorCtx *integratorCtx = ctx;
  struct BLEnsemble *ensemble = integratorCtx->ensemble;
  int numPtcls, i;
  const double complex field = *((const double complex*)x);

  /* For all particles:
   *   compute polarization
   *   compute dipole interaction
   * MPI_Allreduce polarization
   *
   * Scalability can be improved by first computing the polarization,
   * doing an asynchronous all-reduce, then doing the computation of the
   * dipole interaction, and finally finishing the all-reduce.  However,
   * this entails traversing the state arrays twice and evaluating the
   * mode function twice.
   * */
  numPtcls = blRingBufferSize(ensemble->buffer);
  double complex polarization = 0;
  for (i = 0; i < numPtcls; ++i) {
    double mode[3];
    mode[0] = integratorCtx->ex[i];
    mode[1] = integratorCtx->ey[i];
    mode[2] = integratorCtx->ez[i];
    const double complex *psiX =
      (const double complex *)&x[2 + i * ensemble->internalStateSize];
    double complex *psiY =
      (double complex *)&y[2 + i * ensemble->internalStateSize];
    polarization -= I * mode[1] * 1.0e7 * conj(psiX[0]) * psiX[1];
    /* dpsi/dt = -i H psi 
     * H \propto a */
    psiY[0] = -I * mode[1] * 1.0e2 * conj(field) * psiX[1];
    psiY[1] = -I * mode[1] * 1.0e2 * field * psiX[0];
  }
#ifdef BL_WITH_MPI
  MPI_Allreduce(&polarization, y, 2, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#else
  *((double complex*)y) = polarization;
#endif
}

MPI_Request scatterFieldBegin(const struct FieldState *fieldState,
    double *fieldDest) {
#ifdef BL_WITH_MPI
  fieldDest[0] = fieldState->q;
  fieldDest[1] = fieldState->p;
  return 0;
#else
  fieldDest[0] = fieldState->q;
  fieldDest[1] = fieldState->p;
  return 0;
#endif
}

void scatterFieldEnd(MPI_Request req, const struct FieldState *fieldState,
    double *fieldDest) {
#ifdef BL_WITH_MPI
  BL_UNUSED(req);
  BL_UNUSED(fieldState);
  BL_UNUSED(fieldDest);
#else
  BL_UNUSED(req);
  BL_UNUSED(fieldState);
  BL_UNUSED(fieldDest);
#endif
}

static void modeFunction(double x, double y, double z,
                         double *fx, double *fy, double *fz) {
  const double sigmaE = 3.0e-5;
  const double waveNumber = 2.0 * M_PI / 1.0e-6;
  *fx = 0;
  *fy = exp(-(y * y + z * z) / (sigmaE * sigmaE)) * sin(waveNumber * x);
  *fz = 0;
}
