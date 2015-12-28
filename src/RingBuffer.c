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
  int quotient = next / buffer.capacity;
  next -= quotient * buffer.capacity;
  return next;
}

int blRingBufferAddress(struct BLRingBuffer buffer, int index) {
  int address = buffer.begin + index;
  int quotient = address / buffer.capacity;
  address -= quotient * buffer.capacity;
  return address;
}

int blRingBufferPrev(struct BLRingBuffer buffer, int index) {
  int prev = index - 1;
  if (prev < 0) {
    prev += buffer.capacity;
  }
  return prev;
}

int blRingBufferAppendOne(struct BLRingBuffer *buffer) {
  int i = buffer->end;
  buffer->end = blRingBufferNext(*buffer, buffer->end);
  return i;
}

