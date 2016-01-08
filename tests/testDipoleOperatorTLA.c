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
BeforeEach(DipoleOperatorTLA) {
  int i;
  for (i = 0; i < DOF_PER_PTCL * MAX_NUM_PTCLS; ++i) {
    psi[i] = i;
  }
  for (i = 0; i < DOF_PER_PTCL * MAX_NUM_PTCLS; ++i) {
    result[i] = 2.3 * i;
  }
}
AfterEach(DipoleOperatorTLA) {}

Ensure(DipoleOperatorTLA, canBeCreated) {
  dipoleOperator = blDipoleOperatorTLACreate();
  assert_that(dipoleOperator, is_not_null);
  blDipoleOperatorDestroy(dipoleOperator);
}

Ensure(DipoleOperatorTLA, doesNotHaveMatrixElementAlongZ) {
  dipoleOperator = blDipoleOperatorTLACreate();
  ex[0] = 0.0;
  ey[0] = 0.0;
  ez[0] = 1.0;
  blDipoleOperatorApply(dipoleOperator,
                        4, 1, ex, ey, ez, 
                        psi, result, polarization);
  assert_that_double(polarization[0], is_equal_to_double(0.0));
  assert_that_double(polarization[1], is_equal_to_double(0.0));
  int i;
  for (i = 0; i < DOF_PER_PTCL; ++i) {
    assert_that_double(result[i], is_equal_to_double(0.0));
  }
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

