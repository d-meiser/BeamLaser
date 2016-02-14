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
#include <Partition.h>
#include <stdio.h>

static int find_if(int begin, int end, struct BLPredicateClosure pred) {
  while (begin != end && !pred.f(begin, pred.ctx)) {
    ++begin;
  }
  return begin;
}

static int find_backward_if_not(int begin, int end, struct BLPredicateClosure pred) {
  do {
    if (begin == end) {
      return end;
    }
    --end;
  } while (pred.f(end, pred.ctx));
  return ++end;
}

int blPartition(int begin, int end, struct BLPredicateClosure pred,
    struct BLSwapClosure swap) {
  while (1) {
    begin = find_if(begin, end, pred);
    end = find_backward_if_not(begin, end, pred);

    if (begin == end) return begin;

    --end;
    swap.f(begin, end, swap.ctx);
  }
}

static int find_if_bsp(int begin, int end, double pivot, const double *positions) {
  while (begin != end && positions[begin] > pivot) {
    ++begin;
  }
  return begin;
}

static int find_backward_if_not_bsp(int begin, int end, double pivot, const double *positions) {
  do {
    if (begin == end) {
      return end;
    }
    --end;
  } while (positions[end] < pivot);
  return ++end;
}

int blBSP(int begin, int end, double pivot, const double* positions,
    struct BLSwapClosure swap) {
  while (1) {
    begin = find_if_bsp(begin, end, pivot, positions);
    end = find_backward_if_not_bsp(begin, end, pivot, positions);

    if (begin == end) return begin;

    --end;
    swap.f(begin, end, swap.ctx);
  }
}

