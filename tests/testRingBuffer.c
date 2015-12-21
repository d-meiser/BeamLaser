#include <cgreen/cgreen.h>
#include <RingBuffer.h>

Describe(RingBuffer);
BeforeEach(RingBuffer) {}
AfterEach(RingBuffer) {}

Ensure(RingBuffer, canBeConstructedFromSize) {
  struct BLRingBuffer buffer;
  blRingBufferInitialize(3, &buffer);
  assert_that(buffer.capacity, is_equal_to(3));
}

Ensure(RingBuffer, beginsAtZeroAfterConstruction) {
  struct BLRingBuffer buffer;
  blRingBufferInitialize(3, &buffer);
  assert_that(buffer.begin, is_equal_to(0));
}

Ensure(RingBuffer, endsAtZeroAfterConstruction) {
  struct BLRingBuffer buffer;
  blRingBufferInitialize(3, &buffer);
  assert_that(buffer.end, is_equal_to(0));
}

Ensure(RingBuffer, initiallyHasSizeZero) {
  struct BLRingBuffer buffer;
  blRingBufferInitialize(3, &buffer);
  assert_that(blRingBufferSize(buffer), is_equal_to(0));
}

Ensure(RingBuffer, hasCorrectSize) {
  struct BLRingBuffer buffer = {2, 5, 10};
  assert_that(blRingBufferSize(buffer), is_equal_to(3));
}

Ensure(RingBuffer, hasCorrectSizeWithWrapping) {
  struct BLRingBuffer buffer = {5, 1, 10};
  assert_that(blRingBufferSize(buffer), is_equal_to(6));
}

Ensure(RingBuffer, nextIndexIsGreater) {
  struct BLRingBuffer buffer = {5, 1, 10};
  int index = buffer.begin + 2;
  int nextIndex = blRingBufferNext(buffer, index);
  assert_that(nextIndex, is_greater_than(index));
}

Ensure(RingBuffer, nextIndexWrapsCorrectly) {
  struct BLRingBuffer buffer = {5, 1, 10};
  int index = 10;
  int nextIndex = blRingBufferNext(buffer, index);
  assert_that(nextIndex, is_equal_to(0));
}

Ensure(RingBuffer, prevIndexIsSmaller) {
  struct BLRingBuffer buffer = {0, 8, 10};
  int index = buffer.begin + 2;
  int prevIndex = blRingBufferPrev(buffer, index);
  assert_that(prevIndex, is_less_than(index));
}

Ensure(RingBuffer, prevIndexWrapsCorrectly) {
  struct BLRingBuffer buffer = {5, 1, 10};
  int index = 0;
  int prevIndex = blRingBufferPrev(buffer, index);
  assert_that(prevIndex, is_equal_to(9));
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, RingBuffer, canBeConstructedFromSize);
  add_test_with_context(suite, RingBuffer, beginsAtZeroAfterConstruction);
  add_test_with_context(suite, RingBuffer, endsAtZeroAfterConstruction);
  add_test_with_context(suite, RingBuffer, initiallyHasSizeZero);
  add_test_with_context(suite, RingBuffer, hasCorrectSize);
  add_test_with_context(suite, RingBuffer, hasCorrectSizeWithWrapping);
  add_test_with_context(suite, RingBuffer, nextIndexIsGreater);
  add_test_with_context(suite, RingBuffer, nextIndexWrapsCorrectly);
  add_test_with_context(suite, RingBuffer, prevIndexIsSmaller);
  add_test_with_context(suite, RingBuffer, prevIndexWrapsCorrectly);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

