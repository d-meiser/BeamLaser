#include <cgreen/cgreen.h>
#include <Diagnostics.h>
#include <Ensemble.h>
#include <SimulationState.h>

#define MAX_NUM_PTCLS 10
#define INTERNAL_STATE_DIM 5

static struct BLSimulationState simulationState;

Describe(Diagnostics)
BeforeEach(Diagnostics) {
  blEnsembleInitialize(MAX_NUM_PTCLS, INTERNAL_STATE_DIM,
      &simulationState.ensemble);
  int i, j;
  for (i = 0; i < MAX_NUM_PTCLS; ++i) {
    simulationState.ensemble.x[i] = i;
    simulationState.ensemble.y[i] = 2 * i;
    simulationState.ensemble.z[i] = -1.3 * i;
    simulationState.ensemble.vx[i] = 1.0 * i;
    simulationState.ensemble.vy[i] = 3.0 * i;
    simulationState.ensemble.vz[i] = i * i;
    for (j = 0; j < INTERNAL_STATE_DIM; ++j) {
      simulationState.ensemble.internalState[i * INTERNAL_STATE_DIM + j] = i + j;
    }
  }
  simulationState.ensemble.numPtcls = MAX_NUM_PTCLS;
  simulationState.fieldState.q = 1.3;
  simulationState.fieldState.p = 1.7;
}
AfterEach(Diagnostics) {
  blEnsembleFree(&simulationState.ensemble);
}

Ensure(Diagnostics, canBeCreated) {
  struct BLDiagnostics *diagnostics = blDiagnosticsFieldStateCreate(1, 0);
  blDiagnosticsDestroy(diagnostics);
}

Ensure(Diagnostics, worksForMultipleOfDumpPeriodicity) {
  struct BLDiagnostics *diagnostics = blDiagnosticsFieldStateCreate(1, 0);
  blDiagnosticsProcess(diagnostics, 2, &simulationState);
  blDiagnosticsDestroy(diagnostics);
}

Ensure(Diagnostics, doesntDumpWhenNotMultipleOfDumpPeriodicity) {
  struct BLDiagnostics *diagnostics = blDiagnosticsFieldStateCreate(7, 0);
  blDiagnosticsProcess(diagnostics, 5, &simulationState);
  blDiagnosticsDestroy(diagnostics);
}

Ensure(Diagnostics, canCreateMultipleDiagnostics) {
  struct BLDiagnostics *diagnostics = blDiagnosticsFieldStateCreate(7, 0);
  diagnostics = blDiagnosticsPtclsCreate(3, "ptcls", diagnostics);
  diagnostics = blDiagnosticsInternalStateCreate(4, "internal_state", diagnostics);
  blDiagnosticsProcess(diagnostics, 5, &simulationState);
  blDiagnosticsDestroy(diagnostics);
}

Ensure(Diagnostics, canProcessSeveral) {
  struct BLDiagnostics *diagnostics = blDiagnosticsFieldStateCreate(7, 0);
  diagnostics = blDiagnosticsPtclsCreate(7, "ptcls", diagnostics);
  diagnostics = blDiagnosticsInternalStateCreate(7, "internal_state", diagnostics);
  blDiagnosticsProcess(diagnostics, 7, &simulationState);
  blDiagnosticsDestroy(diagnostics);
}

int main(int argn, char **argv)
{
#ifdef BL_WITH_MPI
  MPI_Init(&argn, &argv);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, Diagnostics, canBeCreated);
  add_test_with_context(suite, Diagnostics, worksForMultipleOfDumpPeriodicity);
  add_test_with_context(suite, Diagnostics, doesntDumpWhenNotMultipleOfDumpPeriodicity);
  add_test_with_context(suite, Diagnostics, canCreateMultipleDiagnostics);
  add_test_with_context(suite, Diagnostics, canProcessSeveral);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return result;
}

