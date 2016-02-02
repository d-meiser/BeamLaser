#include <cgreen/cgreen.h>
#include <AtomFieldInteraction.h>


#define MAX_NUM_PTCLS 10
#define DOF_PER_PTCL 2

Describe(AtomFieldInteraction)
BeforeEach(AtomFieldInteraction) {
}
AfterEach(AtomFieldInteraction) {}

Ensure(AtomFieldInteraction, canBeCreated) {
  struct BLAtomFieldInteraction *atomFieldInteraction =
    blAtomFieldInteractionCreate(10, 2, 0, 0);
  assert_that(atomFieldInteraction, is_not_null);
  blAtomFieldInteractionDestroy(atomFieldInteraction);
}

Ensure(AtomFieldInteraction, producesRabiOscillations) {
/* Dipole operator with unit matrix element */
  struct BLDipoleOperator *dipoleOperator = blDipoleOperatorTLACreate(1.0);
  struct BLModeFunction *modeFunction =
    blModeFunctionUniformCreate(0.0, 1.0, 0.0);

  struct BLFieldState fieldState = {1.0, 0.0};

  struct BLEnsemble ensemble;
  blEnsembleInitialize(MAX_NUM_PTCLS, DOF_PER_PTCL, &ensemble);
  int i, j;
  for (i = 0; i < MAX_NUM_PTCLS; ++i) {
    ensemble.x[i] = 0;
    ensemble.y[i] = 0;
    ensemble.z[i] = 0;
    ensemble.vx[i] = 0;
    ensemble.vy[i] = 0;
    ensemble.vz[i] = 0;
    for (j = 0; j < DOF_PER_PTCL; ++j) {
      ensemble.internalState[i * DOF_PER_PTCL + j] = (j == 1) ? 1.0 : 0.0;
    }
  }

  struct BLAtomFieldInteraction *atomFieldInteraction =
    blAtomFieldInteractionCreate(1, DOF_PER_PTCL, dipoleOperator, modeFunction);
  for (i = 0; i < 2; ++i) {
    printf("%le, %le\n",
        creal(ensemble.internalState[i]),
        cimag(ensemble.internalState[i]));
  }

  blAtomFieldInteractionTakeStep(atomFieldInteraction,
      1.0e-3, &fieldState, &ensemble);

  for (i = 0; i < 2; ++i) {
    printf("%le, %le\n",
        creal(ensemble.internalState[i]),
        cimag(ensemble.internalState[i]));
  }

  blEnsembleFree(&ensemble);
  blAtomFieldInteractionDestroy(atomFieldInteraction);
  blModeFunctionDestroy(modeFunction);
  blDipoleOperatorDestroy(dipoleOperator);
}


int main()
{
  TestSuite *suite = create_test_suite();
  add_test_with_context(suite, AtomFieldInteraction, canBeCreated);
  add_test_with_context(suite, AtomFieldInteraction, producesRabiOscillations);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
  return result;
}

