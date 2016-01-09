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
  struct BlBox simulationDomain;
};

struct IntegratorCtx {
  struct BLEnsemble *ensemble;
  struct BLDipoleOperator *dipoleOperator;
  double *ex;
  double *ey;
  double *ez;
};

void setDefaults(struct Configuration *conf);
void particleSink(const struct Configuration *conf, struct BLEnsemble *ensemble);
void processParticleSources(struct ParticleSource *particleSource,
                            struct BLEnsemble *ensemble);
struct ParticleSource *constructParticleSources(
    const struct Configuration *conf);
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
  struct ParticleSource *particleSource;
  BLIntegrator integrator;
  BL_STATUS stat;
  int i, rank;

#ifdef BL_WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = 0;
#endif

  setDefaults(&conf);
  blIntegratorCreate("RK4", conf.maxNumParticles * INTERNAL_STATE_DIM,
                     &integrator);

  stat = blEnsembleInitialize(conf.maxNumParticles, INTERNAL_STATE_DIM,
      &simulationState.ensemble);
  integratorCtx.ensemble = &simulationState.ensemble;
  integratorCtx.dipoleOperator = blDipoleOperatorTLACreate();
  integratorCtx.ex = malloc(conf.maxNumParticles * sizeof(double));
  integratorCtx.ey = malloc(conf.maxNumParticles * sizeof(double));
  integratorCtx.ez = malloc(conf.maxNumParticles * sizeof(double));

  if (stat != BL_SUCCESS) return stat;
  simulationState.fieldState.q = 1.0;
  simulationState.fieldState.p = 0.0;
  if (stat != BL_SUCCESS) return stat;

  particleSource = constructParticleSources(&conf);

  for (i = 0; i < conf.numSteps; ++i) {
    particleSink(&conf, &simulationState.ensemble);
    processParticleSources(particleSource, &simulationState.ensemble);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blFieldAtomInteraction(conf.dt, &simulationState.fieldState,
        &integratorCtx, integrator);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    if (!rank) {
      printf("%d %le %le\n", i,
             simulationState.fieldState.q, simulationState.fieldState.p);
    }
  }


  blParticleSourceDestroy(particleSource);
  blDipoleOperatorDestroy(integratorCtx.dipoleOperator);
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

void particleSink(const struct Configuration *conf,
                  struct BLEnsemble *ensemble) {
  blEnsembleRemoveBelow(conf->simulationDomain.zmin, ensemble->z, ensemble);
}

void processParticleSources(struct ParticleSource *particleSource,
                           struct BLEnsemble *ensemble) {
  int numCreate = blParticleSourceGetNumParticles(particleSource);
  blEnsembleCreateSpace(numCreate, ensemble);
  blParticleSourceCreateParticles(particleSource,
                                  ensemble->x, ensemble->y, ensemble->z,
                                  ensemble->vx, ensemble->vy, ensemble->vz,
                                  ensemble->internalState);
}

void blFieldUpdate(double dt, double kappa, struct FieldState *fieldState) {
  fieldState->q = fieldState->q * exp(-0.5 * kappa * dt);
  fieldState->p = fieldState->p * exp(-0.5 * kappa * dt);
}

void blFieldAtomInteraction(double dt, struct FieldState *fieldState,
        struct IntegratorCtx *integratorCtx, BLIntegrator integrator) {
  int i;
  struct BLEnsemble *ensemble = integratorCtx->ensemble;
  const int numPtcls = ensemble->numPtcls;
  const int fieldOffset = numPtcls * ensemble->internalStateSize;
  int n = numPtcls * ensemble->internalStateSize + 2;

  /* Pack field into internal state buffer. Note that the internalState array
  needs to have room for at least two additional doubles. After the field has
  been scattered to each rank we integrate its equations of motion redundantly.
  */
  MPI_Request fieldRequest =
    scatterFieldBegin(fieldState, ensemble->internalState + fieldOffset);
  for (i = 0; i < ensemble->numPtcls; ++i) {
    modeFunction(ensemble->x[i], ensemble->y[i], ensemble->z[i],
        &integratorCtx->ex[i], &integratorCtx->ey[i], &integratorCtx->ez[i]);
  }
  scatterFieldEnd(fieldRequest, fieldState,
                  ensemble->internalState + fieldOffset);

  blIntegratorTakeStep(integrator, 0.0, dt, n, interactionRHS,
                       ensemble->internalState, ensemble->internalState,
                       integratorCtx);

  fieldState->q = ensemble->internalState[fieldOffset + 0];
  fieldState->p = ensemble->internalState[fieldOffset + 1];
}

void interactionRHS(double t, int n, const double *x, double *y,
                    void *ctx) {
  BL_UNUSED(t);
  BL_UNUSED(n);
  struct IntegratorCtx *integratorCtx = ctx;
  struct BLEnsemble *ensemble = integratorCtx->ensemble;
  int i;
  const int numPtcls = ensemble->numPtcls;
  const int fieldOffset = numPtcls * ensemble->internalStateSize;
  const double complex fieldAmplitude =
    *((const double complex*)(x + fieldOffset));

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
  double complex polarization = 0;
  blDipoleOperatorApply(integratorCtx->dipoleOperator,
                        ensemble->internalStateSize,
                        numPtcls,
                        integratorCtx->ex, integratorCtx->ey, integratorCtx->ez,
                        x, y, (double*)&polarization);

#ifdef BL_WITH_MPI
  MPI_Request polRedReq;
  MPI_Iallreduce(&polarization,
                 y + fieldOffset, 2, MPI_DOUBLE, MPI_SUM,
                 MPI_COMM_WORLD, &polRedReq);
#else
  *((double complex*)(y + fieldOffset)) = polarization;
#endif

  for (i = 0; i < numPtcls; ++i) {
    y[i] *= fieldAmplitude;
  }

#ifdef BL_WITH_MPI
  MPI_Wait(&polRedReq, MPI_STATUS_IGNORE);
#else
#endif
}

MPI_Request scatterFieldBegin(const struct FieldState *fieldState,
    double *fieldDest) {
#ifdef BL_WITH_MPI
  MPI_Request scatterReq;
  MPI_Iscatter(fieldState, 2, MPI_DOUBLE, fieldDest, 2, MPI_DOUBLE, 0,
      MPI_COMM_WORLD, &scatterReq);
  return scatterReq;
#else
  fieldDest[0] = fieldState->q;
  fieldDest[1] = fieldState->p;
  return 0;
#endif
}

void scatterFieldEnd(MPI_Request req, const struct FieldState *fieldState,
    double *fieldDest) {
#ifdef BL_WITH_MPI
  BL_UNUSED(fieldState);
  BL_UNUSED(fieldDest);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
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

struct ParticleSource *constructParticleSources(
    const struct Configuration *conf) {
  struct ParticleSource *particleSource;
  struct BlBox volume = {
    conf->simulationDomain.xmin,
    conf->simulationDomain.xmax,
    conf->simulationDomain.ymin,
    conf->simulationDomain.ymax,
    conf->simulationDomain.zmax,
    conf->simulationDomain.zmax + conf->dt * conf->vbar};
  int numPtcls = round(
      conf->nbar *
      (conf->dt * conf->vbar) /
      (conf->simulationDomain.zmax - conf->simulationDomain.zmin)
      );
  double vbar[] = {0, 0, -conf->vbar};
  double deltaV[] = {
    conf->alpha * conf->deltaV, conf->alpha * conf->deltaV, conf->deltaV};
  double internalState[] = {0, 0, 1, 0};

  particleSource = blParticleSourceUniformCreate(
      volume, numPtcls, vbar, deltaV, 4, internalState, 0);
  return particleSource;
}
