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
#include <Utilities.h>
#include <config.h>
#ifdef BL_WITH_MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>


int main(int argn, char **argv) {
  int result = 0;
#ifdef BL_WITH_MPI
  int status;
  status = MPI_Init(&argn, &argv);
  assert(status == MPI_SUCCESS);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif

  int myRank;
  int numRanks;
#ifdef BL_WITH_MPI
  status = MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  assert(status == MPI_SUCCESS);
  status = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  assert(status == MPI_SUCCESS);
#else
  myRank = 0;
#endif

  double x = 1.0;
  double y = 0;

  BL_MPI_Request req = blAddAllBegin(&x, &y, 1);
  blAddAllEnd(req, &x, &y, 1);

  assert(fabs(y - numRanks * x) < 1.0e-13);

#ifdef BL_WITH_MPI
  status = MPI_Finalize();
  assert(status == MPI_SUCCESS);
#endif
  return result;
}
