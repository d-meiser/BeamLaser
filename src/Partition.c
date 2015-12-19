#include <Partition.h>
#include <stdio.h>

static int find_if(struct BLRingBuffer buffer, struct BLPredicateClosure pred) {
  while (buffer.begin != buffer.end && !pred.f(buffer.begin, pred.ctx)) {
    buffer.begin = blRingBufferNext(buffer, buffer.begin);
  }
  return buffer.begin;
}

static int find_backward_if_not(struct BLRingBuffer buffer, struct BLPredicateClosure pred) {
  do {
    if (buffer.begin == buffer.end) {
      return buffer.end;
    }
    buffer.end = blRingBufferPrev(buffer, buffer.end);
  } while (pred.f(buffer.end, pred.ctx));
  return blRingBufferNext(buffer, buffer.end);
}

int blPartition(struct BLRingBuffer buffer, struct BLPredicateClosure pred,
    struct BLSwapClosure swap) {
  while (1) {
    buffer.begin = find_if(buffer, pred);
    buffer.end = find_backward_if_not(buffer, pred);

    if (buffer.begin == buffer.end) return buffer.begin;

    buffer.end = blRingBufferPrev(buffer, buffer.end);
    swap.f(buffer.begin, buffer.end, swap.ctx);
  }
}
