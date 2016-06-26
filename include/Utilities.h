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
#ifndef UTILITIES_H
#define UTILITIES_H

#include <BeamLaserConfig.h>
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
int blGeneratePoisson(double nbar);

BL_MPI_Request blBcastBegin(const double *src, double *dest, int n);
void blBcastEnd(BL_MPI_Request req, const double *src, double *dest, int n);

BL_MPI_Request blAddAllBegin(const double *src, double *dest, int n);
void blAddAllEnd(BL_MPI_Request req, const double *src, double *dest, int n);

#endif
