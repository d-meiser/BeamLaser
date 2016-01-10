#ifndef UTILITIES_H
#define UTILITIES_H

#include <config.h>
#ifdef BL_WITH_MPI
#include <mpi.h>
#endif

#define BL_UNUSED(a) (void)(a)

#ifdef BL_WITH_MPI
typedef MPI_Request BL_MPI_Request;
#else
typedef int BL_MPI_Request;
#endif


struct BlBox {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
};

double blGenerateGaussianNoise(double mu, double sigma);

BL_MPI_Request blBcastBegin(const double *src, double *dest, int n);
void blBcastEnd(BL_MPI_Request req, const double *src, double *dest, int n);

BL_MPI_Request blAddAllBegin(const double *src, double *dest, int n);
void blAddAllEnd(BL_MPI_Request req, const double *src, double *dest, int n);

#endif
