#include <BeamLaser.h>
#include <stdio.h>
#include <math.h>

static const int MAX_NUM_PARTICLES = 10000;

struct Configuration {
  int numSteps;
  double particleWeight;
  double nbar;
  double dt;
  double vbar;
  double deltaV;
  double alpha;
  struct BBox simulationDomain;
  struct BBox sourceVolume;
};

void setDefaults(struct Configuration *conf);
void computeSourceVolume(struct Configuration *conf);
void particleSink(const struct Configuration *conf, struct BLEnsemble *ensemble);
void particleSource(const struct Configuration *conf, struct BLEnsemble *ensemble);

int main() {
  struct BLEnsemble ensemble;
  struct Configuration conf;
  BL_STATUS stat;
  int i;

  setDefaults(&conf);
  computeSourceVolume(&conf);
  
  stat = blEnsembleInitialize(MAX_NUM_PARTICLES, 4, &ensemble);
  if (stat != BL_SUCCESS) return stat;

  for (i = 0; i < conf.numSteps; ++i) {
    particleSink(&conf, &ensemble);
    particleSource(&conf, &ensemble);
    blEnsemblePush(conf.dt, &ensemble);
    printf("%5d %5d\n", i, blRingBufferSize(ensemble.buffer));
  }

  blEnsembleFree(&ensemble);

  return BL_SUCCESS;
}

void setDefaults(struct Configuration *conf) {
  conf->numSteps = 10;
  conf->particleWeight = 1.0;
  conf->nbar = 1.0e3;
  conf->dt = 1.0e-6;
  conf->vbar = 3.0e2;
  conf->deltaV = 1.0e1;
  conf->alpha = 1.0e-2;
  conf->simulationDomain.xmin = -1.0e-4;
  conf->simulationDomain.xmax = 1.0e-4;
  conf->simulationDomain.ymin = -1.0e-4;
  conf->simulationDomain.ymax = 1.0e-4;
  conf->simulationDomain.zmin = -1.0e-4;
  conf->simulationDomain.zmax = 1.0e-4;
}

void computeSourceVolume(struct Configuration *conf) {
  conf->sourceVolume.xmin = conf->simulationDomain.xmin;
  conf->sourceVolume.xmax = conf->simulationDomain.xmax;
  conf->sourceVolume.ymin = conf->simulationDomain.ymin;
  conf->sourceVolume.ymax = conf->simulationDomain.ymax;
  conf->sourceVolume.zmin = conf->simulationDomain.zmax;
  conf->sourceVolume.zmax = conf->simulationDomain.zmax + conf->vbar * conf->dt;
}

void particleSink(const struct Configuration *conf,
                  struct BLEnsemble *ensemble) {
  blEnsembleRemoveBelow(conf->simulationDomain.zmin, ensemble->z, ensemble);
}

void particleSource(const struct Configuration *conf,
                    struct BLEnsemble *ensemble) {
  int numCreate = round(
      conf->nbar * 
      (conf->simulationDomain.zmax - conf->simulationDomain.zmin) /
      (conf->dt * conf->vbar));
  printf("numCreate: %d\n", numCreate);
  int i;
  while (numCreate > 0) {
    i = blRingBufferAppendOne(&ensemble->buffer);
    blEnsembleCreateParticle(conf->sourceVolume, conf->vbar, conf->deltaV,
                             conf->alpha, i, ensemble);
    --numCreate;
  }
}

