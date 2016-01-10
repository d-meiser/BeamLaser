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
