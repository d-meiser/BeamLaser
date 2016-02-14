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
#include <Update.h>

static struct BLUpdate *update;

Describe(Update)
BeforeEach(Update) {}
AfterEach(Update) {}


Ensure(Update, canBeCreated) {
  update = blUpdateIdentityCreate();
  assert_that(update, is_not_null);
  blUpdateDestroy(update);
}

Ensure(Update, takesStep) {
  update = blUpdateIdentityCreate();
  blUpdateTakeStep(update, 0, 0.1, 0);
  blUpdateDestroy(update);
}

int main(int argn, char **argv)
{
#ifdef BL_WITH_MPI
  MPI_Init(&argn, &argv);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Update, canBeCreated);
  add_test_with_context(suite, Update, takesStep);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return result;
}

