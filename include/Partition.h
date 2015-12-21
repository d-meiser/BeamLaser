#ifndef PARTITION_H
#define PARTITION_H

#include <RingBuffer.h>

struct BLPredicateClosure {
  int (*f)(int i, const void *ctx);
  const void *ctx;
};

struct BLSwapClosure {
  void (*f)(int i, int j, void *ctx);
  void *ctx;
};

int blPartition(struct BLRingBuffer buffer, struct BLPredicateClosure pred,
    struct BLSwapClosure swap);
int blBSP(struct BLRingBuffer buffer, double pivot, const double* positions,
    struct BLSwapClosure swap);

#endif

