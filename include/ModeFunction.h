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
#ifndef MODE_FUNCTION_H
#define MODE_FUNCTION_H

#include <complex.h>

struct BLModeFunction {
  void (*evaluate)(int n,
    const double * restrict x,
    const double * restrict y,
    const double * restrict z,
    double complex * restrict fx,
    double complex * restrict fy,
    double complex * restrict fz,
    void *ctx);
  void (*destroy)(void *);
  void *ctx;
};

void blModeFunctionEvaluate(struct BLModeFunction *modeFunction,
    int n,
    const double * restrict x,
    const double * restrict y,
    const double * restrict z,
    double complex * restrict fx,
    double complex * restrict fy,
    double complex * restrict fz);
void blModeFunctionDestroy(struct BLModeFunction *modeFunction);

struct BLModeFunction *blModeFunctionSimplifiedGaussianCreate(
    double waist, double lambda, double length);

struct BLModeFunction *blModeFunctionUniformCreate(double fx, double fy, double fz);

#endif

