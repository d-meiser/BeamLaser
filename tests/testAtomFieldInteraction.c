#include <cgreen/cgreen.h>
#include <AtomFieldInteraction.h>

Describe(AtomFieldInteraction)
BeforeEach(AtomFieldInteraction) {}
AfterEach(AtomFieldInteraction) {}

Ensure(AtomFieldInteraction, canBeCreated) {
  struct BLAtomFieldInteraction *atomFieldInteraction =
    blAtomFieldInteractionCreate(10, 2, 0, 0);
  blAtomFieldInteractionDestroy(atomFieldInteraction);
}

int main()
{
  TestSuite *suite = create_test_suite();
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

