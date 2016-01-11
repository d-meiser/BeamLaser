#include <BeamLaser.h>
#include <config.h>

#ifdef BL_WITH_MPI
#include <mpi.h>
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <getopt.h>


static const int INTERNAL_STATE_DIM = 4;

struct FieldState {
  double q;
  double p;
};

struct SimulationState {
  double t;
  struct FieldState fieldState;
  struct BLEnsemble ensemble;
};

struct Configuration {
  int numSteps;
  double dt;
  double nbar;
  int maxNumParticles;
  double particleWeight;
  double dipoleMatrixElement;
  double vbar;
  double deltaV;
  double alpha;
  double kappa;
  struct BlBox simulationDomain;
};

struct IntegratorCtx {
  struct BLEnsemble *ensemble;
  struct BLDipoleOperator *dipoleOperator;
  double *ex;
  double *ey;
  double *ez;
};

void setDefaults(struct Configuration *conf);
void processCommandLineArgs(struct Configuration *conf, int argn, char **argv);
void printUsage();
void adjustNumPtclsForNumRanks(struct Configuration *conf);
void particleSink(const struct Configuration *conf, struct BLEnsemble *ensemble);
void processParticleSources(struct ParticleSource *particleSource,
                            struct BLEnsemble *ensemble);
struct ParticleSource *constructParticleSources(
    const struct Configuration *conf);
void blFieldUpdate(double dt, double kappa, struct FieldState *fieldState);
void blFieldAtomInteraction(double dt, struct FieldState *fieldState,
        struct IntegratorCtx *integratorCtx, BLIntegrator integrator);
void interactionRHS(double t, int n, const double *x, double *y,
                    void *ctx);
static void modeFunction(double x, double y, double z,
                         double *fx, double *fy, double *fz);

int main(int argn, char **argv) {
#ifdef BL_WITH_MPI
  MPI_Init(&argn, &argv);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif
  struct SimulationState simulationState;
  struct Configuration conf;
  struct IntegratorCtx integratorCtx;
  struct ParticleSource *particleSource;
  BLIntegrator integrator;
  BL_STATUS stat;
  int i, rank;

#ifdef BL_WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = 0;
#endif

  setDefaults(&conf);
  processCommandLineArgs(&conf, argn, argv);
  adjustNumPtclsForNumRanks(&conf);

  blIntegratorCreate("RK4", conf.maxNumParticles * INTERNAL_STATE_DIM,
                     &integrator);

  stat = blEnsembleInitialize(conf.maxNumParticles, INTERNAL_STATE_DIM,
      &simulationState.ensemble);
  integratorCtx.ensemble = &simulationState.ensemble;
  integratorCtx.dipoleOperator = blDipoleOperatorTLACreate();
  integratorCtx.ex = malloc(conf.maxNumParticles * sizeof(double));
  integratorCtx.ey = malloc(conf.maxNumParticles * sizeof(double));
  integratorCtx.ez = malloc(conf.maxNumParticles * sizeof(double));

  if (stat != BL_SUCCESS) return stat;
  simulationState.fieldState.q = 1.0;
  simulationState.fieldState.p = 0.0;
  if (stat != BL_SUCCESS) return stat;

  particleSource = constructParticleSources(&conf);

  for (i = 0; i < conf.numSteps; ++i) {
    particleSink(&conf, &simulationState.ensemble);
    processParticleSources(particleSource, &simulationState.ensemble);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blFieldAtomInteraction(conf.dt, &simulationState.fieldState,
        &integratorCtx, integrator);
    blFieldUpdate(0.5 * conf.dt, conf.kappa, &simulationState.fieldState);
    blEnsemblePush(0.5 * conf.dt, &simulationState.ensemble);
    if (!rank) {
      printf("%d %le %le\n", i,
             simulationState.fieldState.q, simulationState.fieldState.p);
    }
  }


  blParticleSourceDestroy(particleSource);
  blDipoleOperatorDestroy(integratorCtx.dipoleOperator);
  free(integratorCtx.ex);
  free(integratorCtx.ey);
  free(integratorCtx.ez);
  blIntegratorDestroy(&integrator);
  blEnsembleFree(&simulationState.ensemble);

#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return BL_SUCCESS;
}

void setDefaults(struct Configuration *conf) {
  conf->numSteps = 10;
  conf->particleWeight = 1.0e6;
  conf->dipoleMatrixElement = 1.0e-5 * 1.0e-29;
  conf->nbar = 1.0e3;
  conf->maxNumParticles = 2000;
  conf->dt = 1.0e-8;
  conf->vbar = 3.0e2;
  conf->deltaV = 1.0e1;
  conf->alpha = 1.0e-2;
  conf->kappa = 2.0 * M_PI * 1.0e6;
  conf->simulationDomain.xmin = -1.0e-4;
  conf->simulationDomain.xmax = 1.0e-4;
  conf->simulationDomain.ymin = -1.0e-4;
  conf->simulationDomain.ymax = 1.0e-4;
  conf->simulationDomain.zmin = -1.0e-4;
  conf->simulationDomain.zmax = 1.0e-4;
}

void processCommandLineArgs(struct Configuration *conf, int argn, char **argv) {
  int c;

  while (1)
    {
      static struct option long_options[] =
        {
          {"numSteps",             required_argument, 0, 'n'},
          {"dt",                   required_argument, 0, 'd'},
          {"nbar",                 required_argument, 0, 'N'},
          {"maxNumPtcls",          required_argument, 0, 'm'},
          {"ptclWeight",           required_argument, 0, 'w'},
          {"dipoleMatrixElement",  required_argument, 0, 'D'},
          {"vbar",                 required_argument, 0, 'v'},
          {"deltaV",               required_argument, 0, 'V'},
          {"alpha",                required_argument, 0, 'a'},
          {"kappa",                required_argument, 0, 'K'},
          {"help",                 required_argument, 0, 'h'},
          {0, 0, 0, 0}
        };
      int option_index = 0;

      c = getopt_long(argn, argv, "n:d:N:m:w:D:v:V:a:K:h",
                      long_options, &option_index);

      if (c == -1)
        break;

      switch (c) {
      case 'n':
        if (sscanf(optarg, "%d", &conf->numSteps) != 1) {
          printf("Unable to parse argument to option -n, --numSteps\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'd':
        if (sscanf(optarg, "%lf", &conf->dt) != 1) {
          printf("Unable to parse argument to option -d, --dt\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'N':
        if (sscanf(optarg, "%lf", &conf->nbar) != 1) {
          printf("Unable to parse argument to option -N, --nbar\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'm':
        if (sscanf(optarg, "%d", &conf->maxNumParticles) != 1) {
          printf("Unable to parse argument to option -m, --maxNumPtcls\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'w':
        if (sscanf(optarg, "%lf", &conf->particleWeight) != 1) {
          printf("Unable to parse argument to option -w, --ptclWeight\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'D':
        if (sscanf(optarg, "%lf", &conf->dipoleMatrixElement) != 1) {
          printf("Unable to parse argument to option -D, --dipoleMatrixElement\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'v':
        if (sscanf(optarg, "%lf", &conf->vbar) != 1) {
          printf("Unable to parse argument to option -v, --vbar\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'V':
        if (sscanf(optarg, "%lf", &conf->deltaV) != 1) {
          printf("Unable to parse argument to option -V, --deltaV\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'a':
        if (sscanf(optarg, "%lf", &conf->alpha) != 1) {
          printf("Unable to parse argument to option -a, --alpha\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'K':
        if (sscanf(optarg, "%lf", &conf->kappa) != 1) {
          printf("Unable to parse argument to option -K, --kappa\n");
          printUsage();
          exit(-1);
        }
        break;
      case 'h':
        printUsage();
        exit(0);
      case '?':
        /* getopt_long already printed an error message. */
        break;
      default:
        abort();
      }
    }

  /* Print any remaining command line arguments (not options). */
  if (optind < argn) {
    printf ("non-option ARGV-elements: ");
    while (optind < argn)
      printf ("%s ", argv[optind++]);
    putchar ('\n');
  }
}

void printUsage() {
  printf("\n"
         "BeamLaserTLA --- Simulation of beam laser with two level atoms\n"
         "\n"
         "\n"
         "Usage: BeamLaserTLA [options]\n"
         "\n"
         "Options:\n"
         "-n, --numSteps:           Number of steps to take.\n"
         "-d, --dt:                 Time step size.\n"
         "-N, --nbar:               Mean number of particles in simulation domain.\n"
         "-m, --maxNumPtcl          Maximum number of particles.\n"
         "-w, --ptclWeight          Number of physical particles represented by each\n"
         "                          simulation particle.\n"
         "-D, --dipoleMatrixElement Dipole matrix element of transition with\n"
         "                          Clebsch-Gordan coefficient of one.\n"
         "-v, --vbar                Mean velocity of atoms.\n"
         "-V, --deltaV              Longitudinal velocity spread.\n"
         "-a, --alpha               Beam divergence.\n"
         "-K, --kappa               Cavity damping rate.\n"
         "-h, --help                Print this message.\n"
         "\n"
         );
}

void adjustNumPtclsForNumRanks(struct Configuration *conf) {
  int numRanks = 1;
#ifdef BL_WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
#endif
  conf->nbar /= numRanks;
  conf->maxNumParticles /= numRanks;
}

void particleSink(const struct Configuration *conf,
                  struct BLEnsemble *ensemble) {
  blEnsembleRemoveBelow(conf->simulationDomain.zmin, ensemble->z, ensemble);
}

void processParticleSources(struct ParticleSource *particleSource,
                           struct BLEnsemble *ensemble) {
  int numCreate = blParticleSourceGetNumParticles(particleSource);
  blEnsembleCreateSpace(numCreate, ensemble);
  blParticleSourceCreateParticles(particleSource,
                                  ensemble->x, ensemble->y, ensemble->z,
                                  ensemble->vx, ensemble->vy, ensemble->vz,
                                  ensemble->internalState);
}

void blFieldUpdate(double dt, double kappa, struct FieldState *fieldState) {
  fieldState->q = fieldState->q * exp(-0.5 * kappa * dt);
  fieldState->p = fieldState->p * exp(-0.5 * kappa * dt);
}

void blFieldAtomInteraction(double dt, struct FieldState *fieldState,
        struct IntegratorCtx *integratorCtx, BLIntegrator integrator) {
  int i;
  struct BLEnsemble *ensemble = integratorCtx->ensemble;
  const int numPtcls = ensemble->numPtcls;
  const int fieldOffset = numPtcls * ensemble->internalStateSize;
  int n = numPtcls * ensemble->internalStateSize + 2;

  /* Pack field into internal state buffer. Note that the internalState array
  needs to have room for at least two additional doubles. After the field has
  been scattered to each rank we integrate its equations of motion redundantly.
  */
  BL_MPI_Request fieldRequest =
    blBcastBegin((const double*)fieldState, ensemble->internalState + fieldOffset, 2);
  for (i = 0; i < ensemble->numPtcls; ++i) {
    modeFunction(ensemble->x[i], ensemble->y[i], ensemble->z[i],
        &integratorCtx->ex[i], &integratorCtx->ey[i], &integratorCtx->ez[i]);
  }
  blBcastEnd(fieldRequest, (const double*)fieldState,
                  ensemble->internalState + fieldOffset, 2);

  blIntegratorTakeStep(integrator, 0.0, dt, n, interactionRHS,
                       ensemble->internalState, ensemble->internalState,
                       integratorCtx);

  fieldState->q = ensemble->internalState[fieldOffset + 0];
  fieldState->p = ensemble->internalState[fieldOffset + 1];
}

void interactionRHS(double t, int n, const double *x, double *y,
                    void *ctx) {
  BL_UNUSED(t);
  BL_UNUSED(n);
  struct IntegratorCtx *integratorCtx = ctx;
  struct BLEnsemble *ensemble = integratorCtx->ensemble;
  int i;
  const int numPtcls = ensemble->numPtcls;
  const int fieldOffset = numPtcls * ensemble->internalStateSize;
  const double complex fieldAmplitude =
    *((const double complex*)(x + fieldOffset));

  /* For all particles:
   *   compute polarization
   *   compute dipole interaction
   * MPI_Allreduce polarization
   *
   * Scalability can be improved by first computing the polarization,
   * doing an asynchronous all-reduce, then doing the computation of the
   * dipole interaction, and finally finishing the all-reduce.  However,
   * this entails traversing the state arrays twice and evaluating the
   * mode function twice.
   * */
  double complex polarization = 0;
  blDipoleOperatorApply(integratorCtx->dipoleOperator,
                        ensemble->internalStateSize,
                        numPtcls,
                        integratorCtx->ex, integratorCtx->ey, integratorCtx->ez,
                        x, y, (double*)&polarization);


  BL_MPI_Request polReq =
    blAddAllBegin((const double*)&polarization, y + fieldOffset, 2);
  for (i = 0; i < numPtcls; ++i) {
    y[i] *= fieldAmplitude;
  }
  blAddAllEnd(polReq, (const double*)&polarization, y + fieldOffset, 2);
}


static void modeFunction(double x, double y, double z,
                         double *fx, double *fy, double *fz) {
  const double sigmaE = 3.0e-5;
  const double waveNumber = 2.0 * M_PI / 1.0e-6;
  *fx = 0;
  *fy = exp(-(y * y + z * z) / (sigmaE * sigmaE)) * sin(waveNumber * x);
  *fz = 0;
}

struct ParticleSource *constructParticleSources(
    const struct Configuration *conf) {
  struct ParticleSource *particleSource;
  struct BlBox volume = {
    conf->simulationDomain.xmin,
    conf->simulationDomain.xmax,
    conf->simulationDomain.ymin,
    conf->simulationDomain.ymax,
    conf->simulationDomain.zmax,
    conf->simulationDomain.zmax + conf->dt * conf->vbar};
  int numPtcls = round(
      conf->nbar *
      (conf->dt * conf->vbar) /
      (conf->simulationDomain.zmax - conf->simulationDomain.zmin)
      );
  double vbar[] = {0, 0, -conf->vbar};
  double deltaV[] = {
    conf->alpha * conf->deltaV, conf->alpha * conf->deltaV, conf->deltaV};
  double internalState[] = {0, 0, 1, 0};

  particleSource = blParticleSourceUniformCreate(
      volume, numPtcls, vbar, deltaV, 4, internalState, 0);
  return particleSource;
}

