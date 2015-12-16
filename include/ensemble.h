#ifndef ENSEMBLE_H
#define ENSEMBLE_H

struct Vec3 {
  double x;
  double y;
  double z;
};

struct ComplexNumber {
  double re;
  double im;
};

struct Ensemble {
  int numPtcls;
  int maxNumPtcls;
  struct Vec3 *positions;
  struct Vec3 *velocities;
  struct ComplexNumber *states;
};

#endif

