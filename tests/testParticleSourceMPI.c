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
#include <ParticleSource.h>
#include <config.h>
#ifdef BL_WITH_MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <math.h>


#define MAX_NUM_PTCLS 10


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

  struct BlBox box = {0, 1, 0, 1, 0, 1};
  double vbar[] = {0, 0, 0};
  double deltaV[] = {2.3, 3.5, 7.5};
  double complex initialState[] = {1, 0};
  static const int internalStateSize = 2;
  struct ParticleSource *particleSource =
    blParticleSourceUniformCreate(box, 1, vbar, deltaV,
        internalStateSize, initialState, 0);

  double x[MAX_NUM_PTCLS], y[MAX_NUM_PTCLS], z[MAX_NUM_PTCLS],
         vx[MAX_NUM_PTCLS], vy[MAX_NUM_PTCLS], vz[MAX_NUM_PTCLS];
  double complex internalState[2 * MAX_NUM_PTCLS];

  blParticleSourceCreateParticles(particleSource,
                                  x, y, z, vx, vy, vz,
                                  internalStateSize, internalState);

#ifdef BL_WITH_MPI
  int numRanks;
  status = MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

  /* Compute global center of gravity */
  double meanX, meanY, meanZ;
  status = MPI_Allreduce(&x, &meanX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  assert(status == MPI_SUCCESS);
  status = MPI_Allreduce(&y, &meanY, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  assert(status == MPI_SUCCESS);
  status = MPI_Allreduce(&z, &meanZ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  assert(status == MPI_SUCCESS);

  meanX /= numRanks;
  meanY /= numRanks;
  meanZ /= numRanks;

  /* Make sure that the center of gravity is also in the simulation volume */
  assert(meanX >= box.xmin);
  assert(meanX <= box.xmax);
  assert(meanY >= box.ymin);
  assert(meanY <= box.ymax);
  assert(meanZ >= box.zmin);
  assert(meanZ <= box.zmax);

  /* Make sure the center of gravity doesn't coincide if this ranks center of
  gravity. This is a weak test to make sure that each rank produces independent
  particle data. */
  if (numRanks > 1) {
    assert(fabs(meanX - x[0]) > 1.0e-14);
    assert(fabs(meanY - y[0]) > 1.0e-14);
    assert(fabs(meanZ - z[0]) > 1.0e-14);
  }
#endif


  blParticleSourceDestroy(particleSource);

#ifdef BL_WITH_MPI
  status = MPI_Finalize();
  BL_UNUSED(status);
  assert(status == MPI_SUCCESS);
#endif
  return result;
}
