#include <cgreen/cgreen.h>

Describe(RingBuffer);
BeforeEach(RingBuffer) {}
AfterEach(RingBuffer) {}

Ensure(RingBuffer, canBeConstructedFromSize) {
  assert_that(0, is_equal_to(0));
}

int main(int argc, char **argv)
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, RingBuffer, canBeConstructedFromSize);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

