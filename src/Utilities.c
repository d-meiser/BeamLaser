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
#include <BeamLaserConfig.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define SIMPLE_SPRNG
#include <sprng.h>


static int blGeneratePoissonKnuth(double lambda);


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

int blGeneratePoisson(double nbar) {
  static const double thresholdKnuth = 25.0;
  if (nbar < thresholdKnuth) {
    return blGeneratePoissonKnuth(nbar);
  } else {
    int sample;
    do {
      sample = round(blGenerateGaussianNoise(nbar, sqrt(nbar)));
    } while (sample < 0);
    return sample;
  }
}

static int blGeneratePoissonKnuth(double lambda) {
  double L = exp(-lambda);
  int k = 0;
  double p = 1.0;
  do {
    ++k;
    p *= sprng();
  } while (p > L);
  return k - 1;
}

BL_MPI_Request blBcastBegin(const double *src, double *dest, int n) {
  memcpy(dest, src, n * sizeof(*src));
#ifdef BL_WITH_MPI
  BL_MPI_Request bcastReq;
  MPI_Ibcast(dest, n, MPI_DOUBLE, 0, MPI_COMM_WORLD, &bcastReq);
  return bcastReq;
#else
  return 0;
#endif
}

void blBcastEnd(BL_MPI_Request req, const double *src, double *dest, int n) {
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

BL_MPI_Request blAddAllBegin(const double *src, double *dest, int n) {
#ifdef BL_WITH_MPI
  MPI_Request req;
  MPI_Iallreduce(src, dest, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req);
  return req;
#else
  memcpy(dest, src, n * sizeof(*src));
  return 0;
#endif
}

void blAddAllEnd(BL_MPI_Request req, const double *src, double *dest, int n) {
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
