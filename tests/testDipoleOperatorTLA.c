#include <cgreen/cgreen.h>
#include <DipoleOperator.h>

static struct DipoleOperator *dipoleOperator;
#define MAX_NUM_PTCLS 10
#define DOF_PER_PTCL 4
static double ex[MAX_NUM_PTCLS];
static double ey[MAX_NUM_PTCLS];
static double ez[MAX_NUM_PTCLS];
static double psi[MAX_NUM_PTCLS * DOF_PER_PTCL];
static double result[MAX_NUM_PTCLS * DOF_PER_PTCL];
static double polarization[2];


Describe(DipoleOperatorTLA)
BeforeEach(DipoleOperatorTLA) {}
AfterEach(DipoleOperatorTLA) {}

Ensure(DipoleOperatorTLA, canBeCreated) {
  dipoleOperator = blDipoleOperatorTLACreate();
  assert_that(dipoleOperator, is_not_null);
  blDipoleOperatorDestroy(dipoleOperator);
}

Ensure(DipoleOperatorTLA, doesNotHaveMatrixElementAlongZ) {
  dipoleOperator = blDipoleOperatorTLACreate();
  blDipoleOperatorApply(dipoleOperator,
                        4, 1, ex, ey, ez, 
                        psi, result, polarization);
  blDipoleOperatorDestroy(dipoleOperator);
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, DipoleOperatorTLA, canBeCreated);
  add_test_with_context(suite, DipoleOperatorTLA,
                        doesNotHaveMatrixElementAlongZ);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

