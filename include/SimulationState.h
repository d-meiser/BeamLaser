#ifndef SIMULATION_STATE_H
#define SIMULATION_STATE_H

#include <Ensemble.h>

struct FieldState {
  double q;
  double p;
};

struct SimulationState {
  double t;
  struct FieldState fieldState;
  struct BLEnsemble ensemble;
};

#endif
