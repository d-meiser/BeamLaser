#include <Partition.h>
#include <stdio.h>

static int find_if(int begin, int end, struct BLPredicateClosure pred) {
  while (begin != end && !pred.f(begin, pred.ctx)) {
    ++begin;
  }
  return begin;
}

static int find_backward_if_not(int begin, int end, struct BLPredicateClosure pred) {
  do {
    if (begin == end) {
      return end;
    }
    --end;
  } while (pred.f(end, pred.ctx));
  return ++end;
}

int blPartition(int begin, int end, struct BLPredicateClosure pred,
    struct BLSwapClosure swap) {
  while (1) {
    begin = find_if(begin, end, pred);
    end = find_backward_if_not(begin, end, pred);

    if (begin == end) return begin;

    --end;
    swap.f(begin, end, swap.ctx);
  }
}

static int find_if_bsp(int begin, int end, double pivot, const double *positions) {
  while (begin != end && positions[begin] > pivot) {
    ++begin;
  }
  return begin;
}

static int find_backward_if_not_bsp(int begin, int end, double pivot, const double *positions) {
  do {
    if (begin == end) {
      return end;
    }
    --end;
  } while (positions[end] < pivot);
  return ++end;
}

int blBSP(int begin, int end, double pivot, const double* positions,
    struct BLSwapClosure swap) {
  while (1) {
    begin = find_if_bsp(begin, end, pivot, positions);
    end = find_backward_if_not_bsp(begin, end, pivot, positions);

    if (begin == end) return begin;

    --end;
    swap.f(begin, end, swap.ctx);
  }
}

