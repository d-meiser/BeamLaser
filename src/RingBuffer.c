#include <RingBuffer.h>

void blRingBufferInitialize(int capacity, struct BLRingBuffer *buffer) {
  buffer->begin = 0;
  buffer->end = 0;
  buffer->capacity = capacity;
}

int blRingBufferSize(struct BLRingBuffer buffer) {
  int size = buffer.end - buffer.begin;
  if (size < 0) {
    size += buffer.capacity;
  }
  return size;
}

int blRingBufferNext(struct BLRingBuffer buffer, int index) {
  int next = index + 1;
  if (next >= buffer.capacity) {
    next -= buffer.capacity;
  }
  return next;
}

