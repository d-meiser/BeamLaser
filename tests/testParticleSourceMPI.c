#include <ParticleSource.h>
#include <config.h>
#ifdef BL_WITH_MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <math.h>

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

  struct BBox box = {0, 1, 0, 1, 0, 1};
  double vbar[] = {0, 0, 0};
  double deltaV[] = {2.3, 3.5, 7.5};
  double internalState[] = {1, 0, 0, 0};
  struct ParticleSource *particleSource =
    blParticleSourceUniformCreate(box, 1, vbar, deltaV, 4, internalState, 0);

  double x[1], y[1], z[1], vx[1], vy[1], vz[1];

  blParticleSourceCreateParticles(particleSource,
                                  x, y, z, vx, vy, vz, internalState);

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
  assert(status == MPI_SUCCESS);
#endif
  return result;
}
