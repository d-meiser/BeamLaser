#ifndef SIMULATION_STATE_H
#define SIMULATION_STATE_H

#include <Ensemble.h>

struct BLFieldState {
  double q;
  double p;
};

struct BLSimulationState {
  double t;
  struct BLFieldState fieldState;
  struct BLEnsemble ensemble;
};

#endif
