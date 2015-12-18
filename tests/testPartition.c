#include <cgreen/cgreen.h>
#include <Partition.h>

static int strPredicate_(int i, const void *ctx) {
  const char *str = ctx;
  return str[i] != 'B';
}

static void charSwap_(int i, int j, void *ctx) {
  char *str = ctx;
  if (i != j) {
    char tmp = str[i];
    str[i] = str[j];
    str[j] = str[i];
  }
}

static struct BLPredicateClosure strPredicate;
static struct BLSwapClosure strSwap;
static struct BLRingBuffer buffer;

Describe(Partition);
BeforeEach(Partition) {
  strPredicate.f = strPredicate_;
  strPredicate.ctx = 0;
  strSwap.f = charSwap_;
  strSwap.ctx = 0;
  buffer.begin = 0;
  buffer.end = 4;
  buffer.capacity = 4;
}
AfterEach(Partition) {}


Ensure(Partition, worksForAllBadElements) {
  char *input = "BBBB";
  strPredicate.ctx = input;
  strSwap.ctx = input;
  int partitionPt = blPartition(buffer, strPredicate, strSwap);
  assert_that(partitionPt, is_equal_to(4));
}

Ensure(Partition, givesZeroForAllGoodElements) {
  char *input = "GGGGGGG";
  strPredicate.ctx = input;
  strSwap.ctx = input;
  buffer.end = 7;
  buffer.capacity = 7;
  int partitionPt = blPartition(buffer, strPredicate, strSwap);
  assert_that(partitionPt, is_equal_to(0));
}

int main(int argc, char **argv)
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Partition, worksForAllBadElements);
  add_test_with_context(suite, Partition, givesZeroForAllGoodElements);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

