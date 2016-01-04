#ifndef PARTITION_H
#define PARTITION_H


struct BLPredicateClosure {
  int (*f)(int i, const void *ctx);
  const void *ctx;
};

struct BLSwapClosure {
  void (*f)(int i, int j, void *ctx);
  void *ctx;
};

int blPartition(int begin, int end, struct BLPredicateClosure pred,
    struct BLSwapClosure swap);
int blBSP(int begin, int end, double pivot, const double* positions,
    struct BLSwapClosure swap);

#endif

