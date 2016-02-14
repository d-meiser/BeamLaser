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
#include <Diagnostics.h>
#include <stdlib.h>
#include <stdio.h>


void blDiagnosticsProcess(struct BLDiagnostics *diagnostics, int i,
    struct BLSimulationState *simulationState) {
  while (diagnostics) {
    diagnostics->process(i, simulationState, diagnostics->ctx);
    diagnostics = diagnostics->next;
  }
}

void blDiagnosticsDestroy(struct BLDiagnostics *diagnostics) {
  if (diagnostics) {
    struct BLDiagnostics *next = diagnostics->next;
    diagnostics->destroy(diagnostics->ctx);
    free(diagnostics);
    blDiagnosticsDestroy(next);
  }
}


/*
 * Field diagnostics
 */
struct BLDiagnosticsFieldStateCtx {
  int dumpPeriodicity;
};

static void blDiagnosticsFieldStateProcess(int i,
    struct BLSimulationState* simulationState, void *c) {
  struct BLDiagnosticsFieldStateCtx *ctx = c;
  int rank = 0;
#ifdef BL_WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0 && ((i % ctx->dumpPeriodicity) == 0)) {
    printf("%9d  %9d  %le  %le\n", i,
        simulationState->ensemble.numPtcls,
        simulationState->fieldState.q,
        simulationState->fieldState.p);
  }
}

static void blDiagnosticsFieldStateDestroy(void *ctx) {
  free(ctx);
}

struct BLDiagnostics* blDiagnosticsFieldStateCreate(int dumpPeriodicity,
    struct BLDiagnostics* next) {
  struct BLDiagnostics *this = malloc(sizeof(*this));
  this->process = blDiagnosticsFieldStateProcess;
  this->destroy = blDiagnosticsFieldStateDestroy;
  struct BLDiagnosticsFieldStateCtx *ctx = malloc(sizeof(*ctx));
  ctx->dumpPeriodicity = dumpPeriodicity;
  this->ctx = ctx;
  this->next = next;
  return this;
}

/*
 * Atom diagnostics
 */
struct BLDiagnosticsPtclsCtx {
  int dumpPeriodicity;
  const char *fileName;
};

static void blDiagnosticsPtclsProcess(int i,
    struct BLSimulationState* simulationState, void *c) {
  struct BLDiagnosticsPtclsCtx *ctx = c;
  if ((i % ctx->dumpPeriodicity) == 0) {
    int rank = 0;
#ifdef BL_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    char fileName[1000];
    sprintf(fileName, "%s_%d_%d.txt", ctx->fileName, rank, i);
    FILE *f = fopen(fileName, "w");
    if (!f) return;
    struct BLEnsemble *ensemble = &simulationState->ensemble;
    int j;
    for (j = 0; j < ensemble->numPtcls; ++j) {
      fprintf(f, "%e ", ensemble->x[j]);
    }
    fprintf(f, "\n");
    for (j = 0; j < ensemble->numPtcls; ++j) {
      fprintf(f, "%e ", ensemble->y[j]);
    }
    fprintf(f, "\n");
    for (j = 0; j < ensemble->numPtcls; ++j) {
      fprintf(f, "%e ", ensemble->z[j]);
    }
    fprintf(f, "\n");
    for (j = 0; j < ensemble->numPtcls; ++j) {
      fprintf(f, "%e ", ensemble->vx[j]);
    }
    fprintf(f, "\n");
    for (j = 0; j < ensemble->numPtcls; ++j) {
      fprintf(f, "%e ", ensemble->vy[j]);
    }
    fprintf(f, "\n");
    for (j = 0; j < ensemble->numPtcls; ++j) {
      fprintf(f, "%e ", ensemble->vz[j]);
    }
    fprintf(f, "\n");
    fclose(f);
  }
}

static void blDiagnosticsPtclsDestroy(void *ctx) {
  free(ctx);
}

struct BLDiagnostics* blDiagnosticsPtclsCreate(int dumpPeriodicity,
    const char *fileName, struct BLDiagnostics* next) {
  struct BLDiagnostics *this = malloc(sizeof(*this));
  this->process = blDiagnosticsPtclsProcess;
  this->destroy = blDiagnosticsPtclsDestroy;
  struct BLDiagnosticsPtclsCtx *ctx = malloc(sizeof(*ctx));
  ctx->dumpPeriodicity = dumpPeriodicity;
  ctx->fileName = fileName;
  this->ctx = ctx;
  this->next = next;
  return this;
}

/*
 * Internal state diagnostics
 */
struct BLDiagnosticsInternalStateCtx {
  int dumpPeriodicity;
  const char *fileName;
};

static void blDiagnosticsInternalStateProcess(int i,
    struct BLSimulationState* simulationState, void *c) {
  struct BLDiagnosticsInternalStateCtx *ctx = c;
  if ((i % ctx->dumpPeriodicity) == 0) {
    int rank = 0;
#ifdef BL_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    char fileName[1000];
    sprintf(fileName, "%s_%d_%d.txt", ctx->fileName, rank, i);
    FILE *f = fopen(fileName, "w");
    if (!f) return;
    struct BLEnsemble *ensemble = &simulationState->ensemble;
    int j, k;
    for (j = 0; j < ensemble->numPtcls; ++j) {
      for (k = 0; k < ensemble->internalStateSize; ++k) {
        fprintf(f, "%e %e ",
            creal(ensemble->internalState[j * ensemble->internalStateSize + k]),
            cimag(ensemble->internalState[j * ensemble->internalStateSize + k])
            );
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
}

static void blDiagnosticsInternalStateDestroy(void *ctx) {
  free(ctx);
}

struct BLDiagnostics* blDiagnosticsInternalStateCreate(int dumpPeriodicity,
    const char *fileName, struct BLDiagnostics* next) {
  struct BLDiagnostics *this = malloc(sizeof(*this));
  this->process = blDiagnosticsInternalStateProcess;
  this->destroy = blDiagnosticsInternalStateDestroy;
  struct BLDiagnosticsInternalStateCtx *ctx = malloc(sizeof(*ctx));
  ctx->dumpPeriodicity = dumpPeriodicity;
  ctx->fileName = fileName;
  this->ctx = ctx;
  this->next = next;
  return this;
}

