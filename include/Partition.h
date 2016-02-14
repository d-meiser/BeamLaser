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
#ifndef PARTITION_H
#define PARTITION_H


struct BLPredicateClosure {
  int (*f)(int i, const void *ctx);
  const void *ctx;
};

struct BLSwapClosure {
  void (*f)(int i, int j, void *ctx);
  void *ctx;
};

int blPartition(int begin, int end, struct BLPredicateClosure pred,
    struct BLSwapClosure swap);
int blBSP(int begin, int end, double pivot, const double* positions,
    struct BLSwapClosure swap);

#endif

