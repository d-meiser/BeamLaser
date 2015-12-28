#include <BeamLaser.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef WITH_MPI
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

void setDefaults(struct Configuration *conf);
void computeSourceVolume(struct Configuration *conf);
void particleSink(const struct Configuration *conf, struct BLEnsemble *ensemble);
void particleSource(const struct Configuration *conf, struct BLEnsemble *ensemble);
void blFieldUpdate(double dt, double kappa, struct FieldState *fieldState);
void blFieldAtomInteraction(double dt, struct FieldState *fieldState,
        struct BLEnsemble *ensemble, BLIntegrator integrator);
void scatterFieldEnd(MPI_Request req, const struct FieldState *fieldState,
    double *fieldDest);
MPI_Request scatterFieldBegin(const struct FieldState *fieldState,
    double *fieldDest);
void computeDipoleInteractions(const double *field, double *internalState);
void interactionRHS(double t, int n, const double *x, double *y,
                    void *ctx);

int main() {
  struct SimulationState simulationState;
  struct Configuration conf;
  BLIntegrator integrator;
  BL_STATUS stat;
  int i;

  setDefaults(&conf);
  computeSourceVolume(&conf);
  blIntegratorCreate("RK4", conf.maxNumParticles * INTERNAL_STATE_DIM,
                     &integrator);
  
  stat = blEnsembleInitialize(conf.maxNumParticles, INTERNAL_STATE_DIM,
      &simulationState.ensemble);
  if (stat != BL_SUCCESS) return stat;
  simulationState.fieldState.q = 0.0;
  simulationState.fieldState.p = 0.0;
  if (stat != BL_SUCCESS) return stat;

  for (i = 0; i < conf.numSteps; ++i) {
    particleSink(&conf, &simulationState.ensemble);
    particleSource(&conf, &simulationState.ensemble);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blFieldAtomInteraction(conf.dt, &simulationState.fieldState,
        &simulationState.ensemble, integrator);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
  }

  blIntegratorDestroy(&integrator);
  blEnsembleFree(&simulationState.ensemble);

  return BL_SUCCESS;
}

void setDefaults(struct Configuration *conf) {
  conf->numSteps = 10;
  conf->particleWeight = 1.0;
  conf->nbar = 1.0e3;
  conf->maxNumParticles = 10000;
  conf->dt = 1.0e-7;
  conf->vbar = 3.0e2;
  conf->deltaV = 1.0e1;
  conf->alpha = 1.0e-2;
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
  int i;
  while (numCreate > 0) {
    i = blRingBufferAppendOne(&ensemble->buffer);
    blEnsembleCreateParticle(conf->sourceVolume, conf->vbar, conf->deltaV,
                             conf->alpha, i, ensemble);
    --numCreate;
  }
}

void blFieldUpdate(double dt, double kappa, struct FieldState *fieldState) {
  fieldState->q = fieldState->q * exp(-0.5 * kappa * dt);
  fieldState->p = fieldState->q * exp(-0.5 * kappa * dt);
}

void blFieldAtomInteraction(double dt, struct FieldState *fieldState,
        struct BLEnsemble *ensemble, BLIntegrator integrator) {
  int i, ip;
  int n = blRingBufferSize(ensemble->buffer) * ensemble->internalStateSize + 2;
  double *x = malloc(n * sizeof(double));

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
  scatterFieldEnd(fieldRequest, fieldState, x);

  blIntegratorTakeStep(integrator, 0.0, dt, n, interactionRHS, x, x, ensemble);

  for (i = 0, ip = ensemble->buffer.begin; ip != ensemble->buffer.end;
      ++i, ip = blRingBufferNext(ensemble->buffer, ip)) {
    memcpy(&ensemble->internalState[ip * ensemble->internalStateSize],
           &x[2 + i * ensemble->internalStateSize],
           ensemble->internalStateSize * sizeof(double));
  }

  free(x);
}

void interactionRHS(double t, int n, const double *x, double *y,
                    void *ctx) {
  BL_UNUSED(t);
  BL_UNUSED(n);
  struct BLEnsemble *ensemble = ctx;
  int numPtcls, i;
  struct FieldState *polarization = (struct FieldState*)y;
  struct FieldState *field = (struct FieldState*)x;

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
  polarization->q = 0;
  polarization->p = 0;
  for (i = 0; i < numPtcls; ++i) {
    polarization->q += x[2 + i * ensemble->internalStateSize];
    polarization->p += x[2 + i * ensemble->internalStateSize + 1];
    y[2 + i * ensemble->internalStateSize] = field->q;
  }
}

MPI_Request scatterFieldBegin(const struct FieldState *fieldState,
    double *fieldDest) {
#ifdef WITH_MPI
#else
  fieldDest[0] = fieldState->q;
  fieldDest[1] = fieldState->p;
  return 0;
#endif
}

void scatterFieldEnd(MPI_Request req, const struct FieldState *fieldState,
    double *fieldDest) {
#ifdef WITH_MPI
#else
  BL_UNUSED(req);
  BL_UNUSED(fieldState);
  BL_UNUSED(fieldDest);
#endif
}

void computeDipoleInteractions(const double *field, double *internalState) {
  BL_UNUSED(field);
  BL_UNUSED(internalState);
}
