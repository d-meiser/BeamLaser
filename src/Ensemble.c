#include <Ensemble.h>
#include <stdlib.h>

BL_STATUS blEnsembleInitialize(int capacity, int internalStateSize,
                               struct BLEnsemble *ensemble) {
  ensemble->buffer.begin = 0;
  ensemble->buffer.end = 0;
  ensemble->buffer.capacity = capacity;

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
