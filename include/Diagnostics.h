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
#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <SimulationState.h>
#include <Ensemble.h>


struct BLDiagnostics {
  void (*process)(int i, struct BLSimulationState *simulationState, void *ctx);
  void (*destroy)(void *ctx);
  void *ctx;
  struct BLDiagnostics *next;
};

void blDiagnosticsProcess(struct BLDiagnostics *diagnostics, int i,
    struct BLSimulationState *simulationState);
void blDiagnosticsDestroy(struct BLDiagnostics *diagnostics);

struct BLDiagnostics* blDiagnosticsFieldStateCreate(int dumpPeriodicity,
    struct BLDiagnostics* next);
struct BLDiagnostics* blDiagnosticsPtclsCreate(int dumpPeriodicity,
    const char *filename, struct BLDiagnostics* next);
struct BLDiagnostics* blDiagnosticsInternalStateCreate(int dumpPeriodicity,
    const char *filename, struct BLDiagnostics* next);

#endif
