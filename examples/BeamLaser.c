#include <BeamLaser.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef WITH_MPI
#define MPI_Request int
#endif

static const int INTERNAL_STATE_DIM = 4;

struct FieldState {
  double q;
  double p;
};

struct SimulationState {
  struct FieldState fieldState;
  struct BLEnsemble ensemble;
};

struct RK4WorkSpace {
  double *k1, *k2, *k3, *k4, *tmp;
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
        struct BLEnsemble *ensemble, struct RK4WorkSpace *work);
static BL_STATUS rk4WorkSpaceCreate(int n, struct RK4WorkSpace *work);
void rk4WorkSpaceDestroy(struct RK4WorkSpace *work);
void scatterFieldEnd(MPI_Request req, const struct FieldState *fieldState,
    double *fieldDest);
MPI_Request scatterFieldBegin(const struct FieldState *fieldState,
    double *fieldDest);
MPI_Request gatherPolarizationBegin(const double *internalState,
      double *polarization);
void gatherPolarizationEnd(MPI_Request polarizationRequest,
    const double *internalState, double *polarization);
void computeDipoleInteractions(const double *field, double *internalState);

int main() {
  struct SimulationState simulationState;
  struct Configuration conf;
  struct RK4WorkSpace work;
  BL_STATUS stat;
  int i;

  setDefaults(&conf);
  computeSourceVolume(&conf);
  
  stat = blEnsembleInitialize(conf.maxNumParticles, INTERNAL_STATE_DIM,
      &simulationState.ensemble);
  if (stat != BL_SUCCESS) return stat;
  simulationState.fieldState.q = 0.0;
  simulationState.fieldState.p = 0.0;
  stat = rk4WorkSpaceCreate(conf.maxNumParticles * INTERNAL_STATE_DIM + 2, &work);
  if (stat != BL_SUCCESS) return stat;

  for (i = 0; i < conf.numSteps; ++i) {
    particleSink(&conf, &simulationState.ensemble);
    particleSource(&conf, &simulationState.ensemble);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blFieldAtomInteraction(conf.dt, &simulationState.fieldState,
        &simulationState.ensemble, &work);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    printf("%5d %5d\n", i, blRingBufferSize(simulationState.ensemble.buffer));
  }

  rk4WorkSpaceDestroy(&work);
  blEnsembleFree(&simulationState.ensemble);

  return BL_SUCCESS;
}

void setDefaults(struct Configuration *conf) {
  conf->numSteps = 15;
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
  printf("numCreate: %d\n", numCreate);
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
        struct BLEnsemble *ensemble, struct RK4WorkSpace *work) {
  int i, ip;
  (void)dt;
  for (i = 0, ip = ensemble->buffer.begin; ip != ensemble->buffer.end;
      ++i, ip = blRingBufferNext(ensemble->buffer, ip)) {
    memcpy(&work->tmp[2 + i * INTERNAL_STATE_DIM],
        &ensemble->internalState[ip * INTERNAL_STATE_DIM],
        INTERNAL_STATE_DIM * sizeof(double));
  }

  MPI_Request fieldRequest = scatterFieldBegin(fieldState, work->tmp);
  MPI_Request polarizationRequest = gatherPolarizationBegin(work->tmp + 2,
      work->k1);
  scatterFieldEnd(fieldRequest, fieldState, work->tmp);
  computeDipoleInteractions(work->tmp, work->k1 + 2);
  gatherPolarizationEnd(polarizationRequest, work->tmp, work->k1);
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
  (void)req;
  (void)fieldState;
  (void)fieldDest;
#endif
}

BL_STATUS rk4WorkSpaceCreate(int n, struct RK4WorkSpace *work) {
  work->k1 = malloc(n * sizeof(*work->k1));
  if (!work->k1) return BL_OUT_OF_MEMORY;
  work->k2 = malloc(n * sizeof(*work->k2));
  if (!work->k2) return BL_OUT_OF_MEMORY;
  work->k3 = malloc(n * sizeof(*work->k3));
  if (!work->k3) return BL_OUT_OF_MEMORY;
  work->k4 = malloc(n * sizeof(*work->k4));
  if (!work->k4) return BL_OUT_OF_MEMORY;
  work->tmp = malloc(n * sizeof(*work->tmp));
  if (!work->tmp) return BL_OUT_OF_MEMORY;
  return BL_SUCCESS;
}

void rk4WorkSpaceDestroy(struct RK4WorkSpace *work) {
  free(work->k1);
  free(work->k2);
  free(work->k3);
  free(work->k4);
  free(work->tmp);
}

MPI_Request gatherPolarizationBegin(const double *internalState,
      double *polarization) {
#ifdef WITH_MPI
#else
  (void)internalState;
  polarization[0] = 0;
  polarization[1] = 0;
  return 0;
#endif
}

void gatherPolarizationEnd(MPI_Request polarizationRequest,
    const double *internalState, double *polarization) {
#ifdef WITH_MPI
#else
  (void)polarizationRequest;
  (void)internalState;
  (void)polarization;
#endif
}

void computeDipoleInteractions(const double *field, double *internalState) {
  (void)field;
  (void)internalState;
}
