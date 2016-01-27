#include <cgreen/cgreen.h>
#include <DipoleOperator.h>
#include <complex.h>

static struct BLDipoleOperator *dipoleOperator;
#define MAX_NUM_PTCLS 10
#define DOF_PER_PTCL 2
static double complex ex[MAX_NUM_PTCLS];
static double complex ey[MAX_NUM_PTCLS];
static double complex ez[MAX_NUM_PTCLS];
static double complex psi[MAX_NUM_PTCLS * DOF_PER_PTCL];
static double complex result[MAX_NUM_PTCLS * DOF_PER_PTCL];


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
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  assert_that(dipoleOperator, is_not_null);
  blDipoleOperatorDestroy(dipoleOperator);
}

Ensure(DipoleOperatorTLA, doesNotHaveMatrixElementAlongZ) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  ex[0] = 0.0;
  ey[0] = 0.0;
  ez[0] = 1.0;
  blDipoleOperatorApply(dipoleOperator,
                        2, 1, ex, ey, ez,
                        psi, result);
  int i;
  for (i = 0; i < DOF_PER_PTCL; ++i) {
    assert_that_double(result[i], is_equal_to_double(0.0));
  }
  blDipoleOperatorDestroy(dipoleOperator);
}

Ensure(DipoleOperatorTLA, hasRightMatrixElementAlongX) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  ex[0] = 1.0;
  ey[0] = 3.0;
  ez[0] = -1.7;

  psi[0] = 1.0;
  psi[1] = 0.0;
  blDipoleOperatorApply(dipoleOperator,
                        2, 1, ex, ey, ez,
                        psi, result);
  assert_that_double(creal(result[0]), is_equal_to_double(0.0));
  assert_that_double(cimag(result[0]), is_equal_to_double(0.0));
  assert_that_double(creal(result[1]), is_equal_to_double(1.0));
  assert_that_double(cimag(result[1]), is_equal_to_double(0.0));

  psi[0] = 0.0;
  psi[1] = 1.0;
  blDipoleOperatorApply(dipoleOperator,
                        2, 1, ex, ey, ez,
                        psi, result);
  assert_that_double(creal(result[0]), is_equal_to_double(1.0));
  assert_that_double(cimag(result[0]), is_equal_to_double(0.0));
  assert_that_double(creal(result[1]), is_equal_to_double(0.0));
  assert_that_double(cimag(result[1]), is_equal_to_double(0.0));
  blDipoleOperatorDestroy(dipoleOperator);
}


int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, DipoleOperatorTLA, canBeCreated);
  add_test_with_context(suite, DipoleOperatorTLA,
                        doesNotHaveMatrixElementAlongZ);
  add_test_with_context(suite, DipoleOperatorTLA, hasRightMatrixElementAlongX);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

