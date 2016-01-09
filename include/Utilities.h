#ifndef UTILITIES_H
#define UTILITIES_H

#define BL_UNUSED(a) (void)(a)

struct BBox {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
};

double generateGaussianNoise(double mu, double sigma);

#endif
