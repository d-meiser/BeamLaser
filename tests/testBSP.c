#include <cgreen/cgreen.h>
#include <Partition.h>

static void bspSwap_(int i, int j, void *ctx) {
  double *positions = ctx;
  if (i != j) {
    double tmp = positions[i];
    positions[i] = positions[j];
    positions[j] = tmp;
  }
}

static struct BLSwapClosure bspSwap;
static struct BLRingBuffer buffer;

Describe(BSP)
BeforeEach(BSP) {
  bspSwap.f = bspSwap_;
  bspSwap.ctx = 0;
  buffer.begin = 0;
  buffer.end = 4;
  buffer.capacity = 5;
}
AfterEach(BSP) {}

Ensure(BSP, leavesArrayUnchangedIfAllElementsToLeftOfPivot) {
  double positions[4] = {0, 1, 2, 3};
  bspSwap.ctx = positions;
  int partitionPt = blBSP(buffer, 4.0, positions, bspSwap);
  assert_that(partitionPt, is_equal_to(4));
  assert_that_double(positions[0], is_equal_to_double(0.0));
  assert_that_double(positions[1], is_equal_to_double(1.0));
  assert_that_double(positions[2], is_equal_to_double(2.0));
  assert_that_double(positions[3], is_equal_to_double(3.0));
}

Ensure(BSP, swapsTwoEntriesOnWrongSidesOfPivot) {
  double positions[2] = {1.0, 0.0};
  buffer.begin = 0;
  buffer.end = 2;
  buffer.capacity = 3;

  bspSwap.ctx = positions;
  int partitionPt = blBSP(buffer, 0.5, positions, bspSwap);
  assert_that(partitionPt, is_equal_to(1));
  assert_that_double(positions[0], is_equal_to_double(0.0));
  assert_that_double(positions[1], is_equal_to_double(1.0));
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, BSP, leavesArrayUnchangedIfAllElementsToLeftOfPivot);
  add_test_with_context(suite, BSP, swapsTwoEntriesOnWrongSidesOfPivot);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

