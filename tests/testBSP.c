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
int begin, end, capacity;

Describe(BSP)
BeforeEach(BSP) {
  bspSwap.f = bspSwap_;
  bspSwap.ctx = 0;
  begin = 0;
  end = 4;
  capacity = 5;
}
AfterEach(BSP) {}

Ensure(BSP, leavesArrayUnchangedIfAllElementsToRightOfPivot) {
  double positions[4] = {0, 1, 2, 3};
  bspSwap.ctx = positions;
  int partitionPt = blBSP(begin, end, -1.0, positions, bspSwap);
  assert_that(partitionPt, is_equal_to(4));
  assert_that_double(positions[0], is_equal_to_double(0.0));
  assert_that_double(positions[1], is_equal_to_double(1.0));
  assert_that_double(positions[2], is_equal_to_double(2.0));
  assert_that_double(positions[3], is_equal_to_double(3.0));
}

Ensure(BSP, swapsTwoEntriesOnWrongSidesOfPivot) {
  double positions[2] = {0.0, 1.0};
  begin = 0;
  end = 2;

  bspSwap.ctx = positions;
  int partitionPt = blBSP(begin, end, 0.5, positions, bspSwap);
  assert_that(partitionPt, is_equal_to(1));
  assert_that_double(positions[0], is_equal_to_double(1.0));
  assert_that_double(positions[1], is_equal_to_double(0.0));
}

int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, BSP, leavesArrayUnchangedIfAllElementsToRightOfPivot);
  add_test_with_context(suite, BSP, swapsTwoEntriesOnWrongSidesOfPivot);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

