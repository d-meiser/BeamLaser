#include <cgreen/cgreen.h>
#include <ParticleSource.h>

int main()
{
  TestSuite *suite = create_test_suite();
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

