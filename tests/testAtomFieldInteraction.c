#include <cgreen/cgreen.h>
#include <AtomFieldInteraction.h>
#include <math.h>


#define MAX_NUM_PTCLS 10
#define DOF_PER_PTCL 2

static struct BLDipoleOperator *dipoleOperator;
static struct BLModeFunction *modeFunction;
static struct BLFieldState fieldState;
static struct BLEnsemble ensemble;
static struct BLAtomFieldInteraction *atomFieldInteraction;
static int i, j;

Describe(AtomFieldInteraction)

BeforeEach(AtomFieldInteraction) {
  dipoleOperator = blDipoleOperatorTLACreate(1.0);
  modeFunction = blModeFunctionUniformCreate(0.0, 1.0, 0.0);

  fieldState.q = 1.0;
  fieldState.p = 0.0;

  blEnsembleInitialize(MAX_NUM_PTCLS, DOF_PER_PTCL, &ensemble);
  ensemble.numPtcls = 1;
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
  atomFieldInteraction =
    blAtomFieldInteractionCreate(MAX_NUM_PTCLS, DOF_PER_PTCL, dipoleOperator, modeFunction);
}

AfterEach(AtomFieldInteraction) {
  blEnsembleFree(&ensemble);
  blAtomFieldInteractionDestroy(atomFieldInteraction);
  blModeFunctionDestroy(modeFunction);
  blDipoleOperatorDestroy(dipoleOperator);
}


Ensure(AtomFieldInteraction, canBeCreated) {
  assert_that(atomFieldInteraction, is_not_null);
}

static double nrm_squared(double complex z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

Ensure(AtomFieldInteraction, isNormConserving) {
  for (i = 0; i < 1000; ++i) {
    blAtomFieldInteractionTakeStep(atomFieldInteraction,
        1.0e-3, &fieldState, &ensemble);
  }
  double nrm = 0;
  for (i = 0; i < 2; ++i) {
    nrm += nrm_squared(ensemble.internalState[i]);
  }
  assert_that_double(nrm, is_equal_to_double(1.0));
}

Ensure(AtomFieldInteraction, producesRabiOscillations) {
  /* Use zero weight particles to produce a constant field */
  ensemble.ptclWeight = 0.0;
  for (i = 0; i < 10000; ++i) {
    blAtomFieldInteractionTakeStep(atomFieldInteraction,
        1.0e-3, &fieldState, &ensemble);
  }
  assert_that_double(fieldState.q, is_equal_to_double(1.0));
  assert_that_double(fieldState.p, is_equal_to_double(0.0));

  /* Have to check phase factors here */
  assert_that_double(creal(ensemble.internalState[0]),
      is_equal_to_double(0.0));
  assert_that_double(cimag(ensemble.internalState[0]),
      is_equal_to_double(-sin(10.0)));
  assert_that_double(creal(ensemble.internalState[1]),
      is_equal_to_double(cos(10.0)));
  assert_that_double(cimag(ensemble.internalState[1]),
      is_equal_to_double(0.0));
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
  add_test_with_context(suite, AtomFieldInteraction, canBeCreated);
  add_test_with_context(suite, AtomFieldInteraction, isNormConserving);
  add_test_with_context(suite, AtomFieldInteraction, producesRabiOscillations);
  int result = run_test_suite(suite, create_text_reporter());
  destroy_test_suite(suite);
#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return result;
}

