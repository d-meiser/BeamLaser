#include <Utilities.h>
#include <config.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define SIMPLE_SPRNG
#include <sprng.h>


double blGenerateGaussianNoise(double mu, double sigma) {
  const double epsilon = DBL_EPSILON;
  const double two_pi = 2.0*3.14159265358979323846;

  static double z0, z1;
  static int generate = 0;
  generate = generate == 0 ? 1 : 0;

  if (!generate)
    return z1 * sigma + mu;

  double u1, u2;
  do {
    u1 = sprng();
    u2 = sprng();
  }
  while ( u1 <= epsilon );

  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
  return z0 * sigma + mu;
}

BL_MPI_Request blScatterBegin(const double *src, double *dest, int n) {
#ifdef BL_WITH_MPI
  BL_MPI_Request scatterReq;
  MPI_Iscatter(src, n, MPI_DOUBLE, dest, n, MPI_DOUBLE, 0,
      MPI_COMM_WORLD, &scatterReq);
  return scatterReq;
#else
  memcpy(dest, src, n * sizeof(*src));
  return 0;
#endif
}

void blScatterEnd(BL_MPI_Request req, const double *src, double *dest, int n) {
#ifdef BL_WITH_MPI
  BL_UNUSED(src);
  BL_UNUSED(dest);
  BL_UNUSED(n);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
#else
  BL_UNUSED(req);
  BL_UNUSED(src);
  BL_UNUSED(dest);
  BL_UNUSED(n);
#endif
}

