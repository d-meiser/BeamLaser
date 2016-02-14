/*
Copyright 2014 Dominic Meiser

This file is part of BeamLaser.

BeamLaser is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

BeamLaser is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with BeamLaser.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <cgreen/cgreen.h>
#include <DipoleOperator.h>
#include <complex.h>
#include <math.h>

static struct BLDipoleOperator *dipoleOperator;
#define MAX_NUM_PTCLS 10
#define DOF_PER_PTCL 2
static double complex ex[MAX_NUM_PTCLS];
static double complex ey[MAX_NUM_PTCLS];
static double complex ez[MAX_NUM_PTCLS];
static double complex dx[MAX_NUM_PTCLS];
static double complex dy[MAX_NUM_PTCLS];
static double complex dz[MAX_NUM_PTCLS];
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

Ensure(DipoleOperatorTLA, hasRightMatrixElementAlongY) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  ex[0] = 3.0;
  ey[0] = 1.0;
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

Ensure(DipoleOperatorTLA, hasNoDipoleAlongX) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  psi[0] = 1.0 / sqrt(2.0);
  psi[1] = 1.0 / sqrt(2.0);
  blDipoleOperatorComputeD(dipoleOperator, 2, 1, psi, dx, dy, dz);
  assert_that_double(dx[0], is_equal_to_double(0.0));
  assert_that_double(dz[0], is_equal_to_double(0.0));
}

Ensure(DipoleOperatorTLA, hasRightDipoleAlongY) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  psi[0] = 1.0 / sqrt(2.0);
  psi[1] = 1.0 / sqrt(2.0);
  blDipoleOperatorComputeD(dipoleOperator, 2, 1, psi, dx, dy, dz);
  assert_that_double(dy[0], is_equal_to_double(0.5));
}

Ensure(DipoleOperatorTLA, spinUpHasZeroDipole) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  psi[0] = 1.0;
  psi[1] = 0.0;
  blDipoleOperatorComputeD(dipoleOperator, 2, 1, psi, dx, dy, dz);
  assert_that_double(dx[0], is_equal_to_double(0.0));
  assert_that_double(dy[0], is_equal_to_double(0.0));
  assert_that_double(dz[0], is_equal_to_double(0.0));
}

Ensure(DipoleOperatorTLA, normalizationDoesntMatter) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  psi[0] = 8.0;
  psi[1] = 8.0;
  blDipoleOperatorComputeD(dipoleOperator, 2, 1, psi, dx, dy, dz);
  assert_that_double(dx[0], is_equal_to_double(0.0));
  assert_that_double(dy[0], is_equal_to_double(0.5));
  assert_that_double(dz[0], is_equal_to_double(0.0));
}


int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, DipoleOperatorTLA, canBeCreated);
  add_test_with_context(suite, DipoleOperatorTLA,
                        doesNotHaveMatrixElementAlongZ);
  add_test_with_context(suite, DipoleOperatorTLA, hasRightMatrixElementAlongY);
  add_test_with_context(suite, DipoleOperatorTLA, hasNoDipoleAlongX);
  add_test_with_context(suite, DipoleOperatorTLA, hasRightDipoleAlongY);
  add_test_with_context(suite, DipoleOperatorTLA, spinUpHasZeroDipole);
  add_test_with_context(suite, DipoleOperatorTLA, normalizationDoesntMatter);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

