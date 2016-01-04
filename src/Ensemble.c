#include <Ensemble.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <Partition.h>

#define SIMPLE_SPRNG
#include <sprng.h>

static double generateGaussianNoise(double mu, double sigma);

BL_STATUS blEnsembleInitialize(int capacity, int internalStateSize,
                               struct BLEnsemble *ensemble) {
  ensemble->numPtcls = 0;
  ensemble->maxNumPtcls = capacity;

  ensemble->internalStateSize = internalStateSize;

  ensemble->x = malloc(capacity * sizeof(ensemble->x));
  if (!ensemble->x) { return BL_OUT_OF_MEMORY; }
  ensemble->y = malloc(capacity * sizeof(ensemble->y));
  if (!ensemble->y) { return BL_OUT_OF_MEMORY; }
  ensemble->z = malloc(capacity * sizeof(ensemble->z));
  if (!ensemble->z) { return BL_OUT_OF_MEMORY; }
  ensemble->vx = malloc(capacity * sizeof(ensemble->vx));
  if (!ensemble->vx) { return BL_OUT_OF_MEMORY; }
  ensemble->vy = malloc(capacity * sizeof(ensemble->vy));
  if (!ensemble->vy) { return BL_OUT_OF_MEMORY; }
  ensemble->vz = malloc(capacity * sizeof(ensemble->vz));
  if (!ensemble->vz) { return BL_OUT_OF_MEMORY; }
  ensemble->internalState = malloc(capacity * internalStateSize *
                                   sizeof(ensemble->internalState));
  if (!ensemble->internalState) { return BL_OUT_OF_MEMORY; }
  return BL_SUCCESS;
}

void blEnsembleFree(struct BLEnsemble *ensemble) {
  ensemble->numPtcls = 0;
  ensemble->maxNumPtcls = 0;

  free(ensemble->x);
  free(ensemble->y);
  free(ensemble->z);
  free(ensemble->vx);
  free(ensemble->vy);
  free(ensemble->vz);
  free(ensemble->internalState);
}

#define SWAP(a,b) do { \
  tmp = (a); \
  (a) = (b); \
  (b) = tmp; \
} while (0)

static void ensembleSwap_(int i, int j, void *ctx) {
  int m;
  struct BLEnsemble *ensemble = ctx;
  if (i != j) {
    double tmp;
    SWAP(ensemble->x[i], ensemble->x[j]);
    SWAP(ensemble->y[i], ensemble->y[j]);
    SWAP(ensemble->z[i], ensemble->z[j]);
    SWAP(ensemble->vx[i], ensemble->vx[j]);
    SWAP(ensemble->vy[i], ensemble->vy[j]);
    SWAP(ensemble->vz[i], ensemble->vz[j]);
    for (m = 0; m < ensemble->internalStateSize; ++m) {
      SWAP(ensemble->internalState[ensemble->internalStateSize * i + m],
           ensemble->internalState[ensemble->internalStateSize * j + m]);
    }
  }
}

void blEnsembleRemoveBelow(double cutoff, double *positions,
                           struct BLEnsemble *ensemble) {
  struct BLSwapClosure swap;
  swap.f = ensembleSwap_;
  swap.ctx = ensemble;
  ensemble->numPtcls = blBSP(0, ensemble->numPtcls, cutoff, positions, swap);
}

void blEnsemblePush(double dt, struct BLEnsemble *ensemble) {
  int i;
  for (i = 0; i < ensemble->numPtcls; ++i) {
    ensemble->x[i] += dt * ensemble->vx[i];
    ensemble->y[i] += dt * ensemble->vy[i];
    ensemble->z[i] += dt * ensemble->vz[i];
  }
}

void blEnsembleCreateParticle(struct BBox box, double vbar, double deltaV,
                              double alpha,
                              int i, struct BLEnsemble* ensemble) {
  ensemble->x[i] = box.xmin + (box.xmax - box.xmin) * sprng();
  ensemble->y[i] = box.ymin + (box.ymax - box.ymin) * sprng();
  ensemble->z[i] = box.zmin + (box.zmax - box.zmin) * sprng();
  ensemble->vx[i] = generateGaussianNoise(0.0, alpha * deltaV);
  ensemble->vy[i] = generateGaussianNoise(0.0, alpha * deltaV);
  ensemble->vz[i] = generateGaussianNoise(-vbar, alpha * deltaV);
}

static double generateGaussianNoise(double mu, double sigma) {
  const double epsilon = DBL_EPSILON;
  const double two_pi = 2.0*3.14159265358979323846;

  static double z0, z1;
  static int generate = 0;
  generate = generate == 0 ? 1 : 0;

  if (!generate)
    return z1 * sigma + mu;

  double u1, u2;
  do {
    u1 = sprng();
    u2 = sprng();
  }
  while ( u1 <= epsilon );

  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
  return z0 * sigma + mu;
}

