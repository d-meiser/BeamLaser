#ifndef SIMULATION_STATE_H
#define SIMULATION_STATE_H

#include <Ensemble.h>
#include <FieldState.h>

struct BLSimulationState {
  double t;
  struct BLFieldState fieldState;
  struct BLEnsemble ensemble;
};

#endif
