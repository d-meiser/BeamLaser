#include <cgreen/cgreen.h>
#include <Partition.h>

static int samplePredicate_(int i, const void *ctx) {
  const int *arr = ctx;
  return arr[i];
}

static void sampleSwap_(int i, int j, void *ctx) {
  int *arr = ctx;
  if (i != j) {
    int tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
  }
}

static struct BLPredicateClosure samplePredicate;
static struct BLSwapClosure sampleSwap;
static struct BLRingBuffer buffer;

Describe(Partition);
BeforeEach(Partition) {
  samplePredicate.f = samplePredicate_;
  samplePredicate.ctx = 0;
  sampleSwap.f = sampleSwap_;
  sampleSwap.ctx = 0;
  buffer.begin = 0;
  buffer.end = 4;
  buffer.capacity = 4;
}
AfterEach(Partition) {}


Ensure(Partition, worksForAllBadElements) {
  int arr[4] = {0, 0, 0, 0};
  samplePredicate.ctx = arr;
  sampleSwap.ctx = arr;
  int partitionPt = blPartition(buffer, samplePredicate, sampleSwap);
  assert_that(partitionPt, is_equal_to(4));
}

Ensure(Partition, givesZeroForAllGoodElements) {
  int arr[4] = {1, 1, 1, 1};
  samplePredicate.ctx = arr;
  sampleSwap.ctx = arr;
  int partitionPt = blPartition(buffer, samplePredicate, sampleSwap);
  assert_that(partitionPt, is_equal_to(0));
}

Ensure(Partition, leavesPartitionedSequenceUnchanged) {
  int arr[7] = {0, 0, 0, 1, 1, 1, 1};
  samplePredicate.ctx = arr;
  sampleSwap.ctx = arr;
  buffer.begin = 0;
  buffer.end = 7;
  buffer.capacity = 7;
  int partitionPt = blPartition(buffer, samplePredicate, sampleSwap);
  assert_that(partitionPt, is_equal_to(3));
}

Ensure(Partition, swapsTwoElementsInWrongOrder) {
  int arr[2] = {1, 0};
  samplePredicate.ctx = arr;
  sampleSwap.ctx = arr;
  buffer.begin = 0;
  buffer.end = 2;
  buffer.capacity = 2;
  int partitionPt = blPartition(buffer, samplePredicate, sampleSwap);
  assert_that(partitionPt, is_equal_to(1));
  assert_that(arr[0], is_equal_to(0));
  assert_that(arr[1], is_equal_to(1));
}

Ensure(Partition, putsBadElementsAtBeginning) {
  /*            0  1  2  3  4  5  6  7 */
  int arr[8] = {1, 0, 0, 1, 1, 1, 0, 1};
  samplePredicate.ctx = arr;
  sampleSwap.ctx = arr;
  buffer.begin = 0;
  buffer.end = 8;
  buffer.capacity = 8;
  int partitionPt = blPartition(buffer, samplePredicate, sampleSwap);
  assert_that(partitionPt, is_equal_to(3));
  assert_that(arr[0], is_equal_to(0));
  assert_that(arr[1], is_equal_to(0));
  assert_that(arr[2], is_equal_to(0));
  assert_that(arr[3], is_equal_to(1));
  assert_that(arr[4], is_equal_to(1));
  assert_that(arr[5], is_equal_to(1));
  assert_that(arr[6], is_equal_to(1));
  assert_that(arr[7], is_equal_to(1));
}

int main(int argc, char **argv)
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Partition, worksForAllBadElements);
  add_test_with_context(suite, Partition, givesZeroForAllGoodElements);
  add_test_with_context(suite, Partition, leavesPartitionedSequenceUnchanged);
  add_test_with_context(suite, Partition, swapsTwoElementsInWrongOrder);
  add_test_with_context(suite, Partition, putsBadElementsAtBeginning);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

