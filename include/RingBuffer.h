#ifndef RING_BUFFER_H
#define RING_BUFFER_H

struct BLRingBuffer {
  int begin;
  int end;
  int capacity;
};

void blRingBufferInitialize(int capacity, struct BLRingBuffer *buffer);
int blRingBufferSize(struct BLRingBuffer buffer);
int blRingBufferNext(struct BLRingBuffer buffer, int index);
int blRingBufferPrev(struct BLRingBuffer buffer, int index);
int blRingBufferAppendOne(struct BLRingBuffer *buffer);

#endif

