#include <Utilities.h>
#include <config.h>
#include <math.h>
#include <float.h>
#define SIMPLE_SPRNG
#include <sprng.h>


double generateGaussianNoise(double mu, double sigma) {
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

