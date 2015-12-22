#include <Ensemble.h>
#include <stdlib.h>
#include <Partition.h>

BL_STATUS blEnsembleInitialize(int capacity, int internalStateSize,
                               struct BLEnsemble *ensemble) {
  ensemble->buffer.begin = 0;
  ensemble->buffer.end = 0;
  ensemble->buffer.capacity = capacity;

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
  ensemble->buffer.begin = 0;
  ensemble->buffer.end = 0;
  ensemble->buffer.capacity = 0;

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
  ensemble->buffer.begin = blBSP(ensemble->buffer, cutoff, positions, swap);
}
