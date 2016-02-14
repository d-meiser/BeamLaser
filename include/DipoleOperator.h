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
#ifndef DIPOLE_OPERATOR_H
#define DIPOLE_OPERATOR_H

#include <complex.h>


struct BLDipoleOperator {
  void (*apply)(int stride, int numPtcls,
      const double complex *ex,
      const double complex *ey,
      const double complex *ez,
      const double complex *psi, double complex *result, void *ctx);
  void (*computeD)(int stride, int numPtcls,
      const double complex *psi,
      double complex *dx, double complex *dy, double complex *dz,
      void *ctx);
  void (*destroy)(void *ctx);
  void *ctx;
};

void blDipoleOperatorApply(struct BLDipoleOperator *op,
    int stride, int numPtcls,
    const double complex *ex,
    const double complex *ey,
    const double complex *ez,
    const double complex *psi, double complex *result);

void blDipoleOperatorComputeD(struct BLDipoleOperator *op,
    int stride, int numPtcls, const double complex *psi,
    double complex *dx, double complex *dy, double complex *dz);

void blDipoleOperatorExpectationValue(struct BLDipoleOperator *op,
    int stride, int numPtcls,
    const double complex *psi,
    double complex *dx,
    double complex *dy,
    double complex *dz);

void blDipoleOperatorDestroy(struct BLDipoleOperator *op);


struct BLDipoleOperator *blDipoleOperatorTLACreate(double dipoleMatrixElement);

#endif

