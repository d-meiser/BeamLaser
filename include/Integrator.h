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
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Errors.h>


typedef struct BLIntegrator_ *BLIntegrator;

typedef void (*BLIntegratorRHS)(double t, int n, const double *x, double *y,
    void *ctx);

BL_STATUS blIntegratorCreate(const char* name, int n, BLIntegrator *integrator);
void blIntegratorDestroy(BLIntegrator *integrator);
void blIntegratorTakeStep(BLIntegrator integrator, double t, double dt, int n,
    BLIntegratorRHS rhs, const double *x, double *y, void *ctx);

#endif
