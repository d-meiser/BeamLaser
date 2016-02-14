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


static const int INTERNAL_STATE_DIM = 2;

struct Configuration {
  int numSteps;
  int dumpPeriod;
  int dumpField;
  int dumpPhaseSpace;
  int dumpInternalState;
  double dt;
  double nbar;
  int maxNumParticles;
  double dipoleMatrixElement;
  double ptclWeight;
  double vbar;
  double deltaV;
  double alpha;
  double kappa;
  int uniform;
  double waist;
  double lambda;
  double length;
  struct BlBox simulationDomain;
};


void setDefaults(struct Configuration *conf);
void processCommandLineArgs(struct Configuration *conf, int argn, char **argv);
void printUsage(const char* errorMessage);
void adjustNumPtclsForNumRanks(struct Configuration *conf);
void processParticleSources(struct BLParticleSource *particleSource,
                            struct BLEnsemble *ensemble);
struct BLParticleSource *constructParticleSources(
    const struct Configuration *conf);
struct BLDiagnostics *constructDiagnostics(const struct Configuration *conf);
struct BLModeFunction *constructModeFunction( const struct Configuration *conf);
void interactionRHS(double t, int n, const double *x, double *y,
                    void *ctx);


int main(int argn, char **argv) {
#ifdef BL_WITH_MPI
  MPI_Init(&argn, &argv);
#else
  BL_UNUSED(argn);
  BL_UNUSED(argv);
#endif
  struct BLSimulationState simulationState;
  struct Configuration conf;
  struct BLParticleSource *particleSource;
  struct BLDipoleOperator *dipoleOperator;
  struct BLModeFunction *modeFunction;
  struct BLDiagnostics *diagnostics = 0;
  BL_STATUS stat;
  int i, rank;

#ifdef BL_WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = 0;
  BL_UNUSED(rank);
#endif

  setDefaults(&conf);
  processCommandLineArgs(&conf, argn, argv);
  adjustNumPtclsForNumRanks(&conf);

  stat = blEnsembleCreate(conf.maxNumParticles, INTERNAL_STATE_DIM,
      &simulationState.ensemble);
  simulationState.ensemble.ptclWeight = conf.ptclWeight;

  if (stat != BL_SUCCESS) return stat;
  simulationState.fieldState.q = 1.0;
  simulationState.fieldState.p = 0.0;
  if (stat != BL_SUCCESS) return stat;

  particleSource = constructParticleSources(&conf);
  diagnostics = constructDiagnostics(&conf);
  dipoleOperator = blDipoleOperatorTLACreate(conf.dipoleMatrixElement);
  modeFunction = constructModeFunction(&conf);

  struct BLUpdate *fieldUpdate = blFieldUpdateCreate(conf.kappa, conf.kappa);
  struct BLUpdate *atomPush = blPushUpdateCreate();
  struct BLUpdate *atomFieldInteraction = blAtomFieldInteractionCreate(
      simulationState.ensemble.maxNumPtcls,
      simulationState.ensemble.internalStateSize,
      dipoleOperator, modeFunction);
  struct BLUpdate *sinks = blSinkBelowCreate(conf.simulationDomain.zmin);

  for (i = 0; i < conf.numSteps; ++i) {
    blUpdateTakeStep(sinks, i * conf.dt, conf.dt, &simulationState);
    processParticleSources(particleSource, &simulationState.ensemble);
    blUpdateTakeStep(atomPush, i * conf.dt, 0.5 * conf.dt, &simulationState);
    blUpdateTakeStep(fieldUpdate, i * conf.dt, 0.5 * conf.dt, &simulationState);
    blUpdateTakeStep(atomFieldInteraction, i * conf.dt, conf.dt, &simulationState);
    blUpdateTakeStep(fieldUpdate, (i + 1) * conf.dt, 0.5 * conf.dt, &simulationState);
    blUpdateTakeStep(atomPush, (i + 1) * conf.dt, 0.5 * conf.dt, &simulationState);
    blDiagnosticsProcess(diagnostics, i, &simulationState);
  }

  blUpdateDestroy(sinks);
  blUpdateDestroy(atomFieldInteraction);
  blUpdateDestroy(atomPush);
  blUpdateDestroy(fieldUpdate);
  blModeFunctionDestroy(modeFunction);
  blDipoleOperatorDestroy(dipoleOperator);
  blDiagnosticsDestroy(diagnostics);
  blParticleSourceDestroy(particleSource);
  blEnsembleDestroy(&simulationState.ensemble);

#ifdef BL_WITH_MPI
  MPI_Finalize();
#endif
  return BL_SUCCESS;
}

void setDefaults(struct Configuration *conf) {
  conf->numSteps = 10;
  conf->dumpPeriod = 1;
  conf->dumpField = 0;
  conf->dumpPhaseSpace = 0;
  conf->dumpInternalState = 0;
  conf->dipoleMatrixElement = 1.0e-29;
  conf->ptclWeight = 1.0e0;
  conf->nbar = 1.0e3;
  conf->maxNumParticles = 2000;
  conf->dt = 1.0e-8;
  conf->vbar = 3.0e2;
  conf->deltaV = 1.0e1;
  conf->alpha = 1.0e-2;
  conf->kappa = 2.0 * M_PI * 1.0e6;
  conf->simulationDomain.xmin = -5.0e-5;
  conf->simulationDomain.xmax = 5.0e-5;
  conf->simulationDomain.ymin = -5.0e-5;
  conf->simulationDomain.ymax = 5.0e-5;
  conf->simulationDomain.zmin = -1.0e-4;
  conf->simulationDomain.zmax = 1.0e-4;
  conf->uniform = 0;
  conf->waist = 1.0e-4;
  conf->lambda = 1.0e-6;
  conf->length = 1.0e-2;
}

void processCommandLineArgs(struct Configuration *conf, int argn, char **argv) {
  int c;
  double tmp;

  while (1)
    {
      static struct option long_options[] =
        {
          {"numSteps",             required_argument, 0, 'n'},
          {"dumpPeriod",           required_argument, 0, 'p'},
          {"dumpField",            no_argument,       0, 'f'},
          {"dumpPhaseSpace",       no_argument,       0, 's'},
          {"dumpInternalState",    no_argument,       0, 'i'},
          {"dt",                   required_argument, 0, 'd'},
          {"nbar",                 required_argument, 0, 'N'},
          {"maxNumPtcls",          required_argument, 0, 'm'},
          {"ptclWeight",           required_argument, 0, 'w'},
          {"dipoleMatrixElement",  required_argument, 0, 'D'},
          {"vbar",                 required_argument, 0, 'v'},
          {"deltaV",               required_argument, 0, 'V'},
          {"alpha",                required_argument, 0, 'a'},
          {"kappa",                required_argument, 0, 'K'},
          {"help",                 no_argument,       0, 'h'},
          {"uniform",              no_argument,       0, 'u'},
          {"waist",                required_argument, 0, 'W'},
          {"lambda",               required_argument, 0, 'l'},
          {"length",               required_argument, 0, 'L'},
          {0, 0, 0, 0}
        };
      int option_index = 0;

      c = getopt_long(argn, argv, "n:p:fsid:N:m:w:D:v:V:a:K:huw:l:L:",
                      long_options, &option_index);

      if (c == -1)
        break;

      switch (c) {
      case 'n':
        if (sscanf(optarg, "%lf", &tmp) != 1) {
          printUsage("Unable to parse argument to option -n, --numSteps\n");
          exit(-1);
        }
        conf->numSteps = tmp;
        break;
      case 'p':
        if (sscanf(optarg, "%lf", &tmp) != 1) {
          printUsage("Unable to parse argument to option -p, --dumpPeriod\n");
          exit(-1);
        }
        conf->dumpPeriod = tmp;
        break;
      case 'f':
        conf->dumpField = 1;
        break;
      case 's':
        conf->dumpPhaseSpace = 1;
        break;
      case 'i':
        conf->dumpInternalState = 1;
        break;
      case 'd':
        if (sscanf(optarg, "%lf", &conf->dt) != 1) {
          printUsage("Unable to parse argument to option -d, --dt\n");
          exit(-1);
        }
        break;
      case 'N':
        if (sscanf(optarg, "%lf", &conf->nbar) != 1) {
          printUsage("Unable to parse argument to option -N, --nbar\n");
          exit(-1);
        }
        break;
      case 'm':
        if (sscanf(optarg, "%lf", &tmp) != 1) {
          printUsage("Unable to parse argument to option -m, --maxNumPtcls\n");
          exit(-1);
        }
        conf->maxNumParticles = (int)tmp;
        break;
      case 'w':
        if (sscanf(optarg, "%lf", &conf->ptclWeight) != 1) {
          printUsage("Unable to parse argument to option -w, --ptclWeight\n");
          exit(-1);
        }
        break;
      case 'D':
        if (sscanf(optarg, "%lf", &conf->dipoleMatrixElement) != 1) {
          printUsage("Unable to parse argument to option -D, --dipoleMatrixElement\n");
          exit(-1);
        }
        break;
      case 'v':
        if (sscanf(optarg, "%lf", &conf->vbar) != 1) {
          printUsage("Unable to parse argument to option -v, --vbar\n");
          exit(-1);
        }
        break;
      case 'V':
        if (sscanf(optarg, "%lf", &conf->deltaV) != 1) {
          printUsage("Unable to parse argument to option -V, --deltaV\n");
          exit(-1);
        }
        break;
      case 'a':
        if (sscanf(optarg, "%lf", &conf->alpha) != 1) {
          printUsage("Unable to parse argument to option -a, --alpha\n");
          exit(-1);
        }
        break;
      case 'K':
        if (sscanf(optarg, "%lf", &conf->kappa) != 1) {
          printUsage("Unable to parse argument to option -K, --kappa\n");
          exit(-1);
        }
        break;
      case 'u':
        conf->uniform = 1;
        break;
      case 'W':
        if (sscanf(optarg, "%lf", &conf->waist) != 1) {
          printUsage("Unable to parse argument to option -W, --waist\n");
          exit(-1);
        }
        break;
      case 'l':
        if (sscanf(optarg, "%lf", &conf->lambda) != 1) {
          printUsage("Unable to parse argument to option -l, --lambda\n");
          exit(-1);
        }
        break;
      case 'L':
        if (sscanf(optarg, "%lf", &conf->length) != 1) {
          printUsage("Unable to parse argument to option -L, --length\n");
          exit(-1);
        }
        break;
      case 'h':
        printUsage(0);
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

void printUsage(const char* errorMessage) {
  int myRank = 0;
#ifdef BL_WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif
  if (myRank) return;
  if (errorMessage) {
    printf("%s\n", errorMessage);
  }
  printf("\n"
           "BeamLaserTLA --- Simulation of beam laser with two level atoms\n"
           "\n"
           "\n"
           "Usage: BeamLaserTLA [options]\n"
           "\n"
           "Options:\n"
           "-n, --numSteps:           Number of steps to take.\n"
           "-p, --dumpPeriod:         Number of steps between dumps.\n"
           "-f, --dumpField:          Whether to dump the field.\n"
           "-s, --dumpPhaseSpace:     Whether to dump phase space coordinates of the particles.\n"
           "-i, --dumpInternalState:  Whether to dump the particles' internal state.\n"
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
           "-u, --uniform             Use uniform mode function.\n"
           "-w, --waist               1/e radius of mode.\n"
           "-l, --lambda              Wave length of the mode.\n"
           "-L, --length              Length of the mode.\n"
           "-h, --help                Print this message.\n"
           "\n"
           );
  exit(0);
}

void adjustNumPtclsForNumRanks(struct Configuration *conf) {
  int numRanks = 1;
#ifdef BL_WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
#endif
  conf->nbar /= numRanks;
  conf->maxNumParticles /= numRanks;
}

void processParticleSources(struct BLParticleSource *particleSource,
                           struct BLEnsemble *ensemble) {
  int numCreate = blParticleSourceGetNumParticles(particleSource);
  blEnsembleCreateSpace(numCreate, ensemble);
  blParticleSourceCreateParticles(particleSource,
                                  ensemble->x, ensemble->y, ensemble->z,
                                  ensemble->vx, ensemble->vy, ensemble->vz,
                                  ensemble->internalStateSize,
                                  ensemble->internalState);
}

struct BLParticleSource *constructParticleSources(
    const struct Configuration *conf) {
  struct BLParticleSource *particleSource;
  struct BlBox volume = {
    conf->simulationDomain.xmin,
    conf->simulationDomain.xmax,
    conf->simulationDomain.ymin,
    conf->simulationDomain.ymax,
    conf->simulationDomain.zmax,
    conf->simulationDomain.zmax + conf->dt * conf->vbar};
  double nbar = conf->nbar *
      (conf->dt * conf->vbar) /
      (conf->simulationDomain.zmax - conf->simulationDomain.zmin);
  double vbar[] = {0, 0, -conf->vbar};
  double deltaV[] = {
    conf->alpha * conf->deltaV, conf->alpha * conf->deltaV, conf->deltaV};
  double complex internalState[2];
  internalState[0] = 0;
  internalState[1] = 1;

  particleSource = blParticleSourceUniformCreate(
      volume, nbar, vbar, deltaV, 2, internalState, 0);
  return particleSource;
}

struct BLDiagnostics *constructDiagnostics(
    const struct Configuration *conf) {
  struct BLDiagnostics *diagnostics = 0;
  if (conf->dumpField) {
    diagnostics = blDiagnosticsFieldStateCreate(conf->dumpPeriod, diagnostics);
  }
  if (conf->dumpPhaseSpace) {
    diagnostics = blDiagnosticsPtclsCreate(conf->dumpPeriod, "ptcls",
        diagnostics);
  }
  if (conf->dumpInternalState) {
    diagnostics = blDiagnosticsInternalStateCreate(conf->dumpPeriod,
        "internal_state", diagnostics);
  }
  return diagnostics;
}

struct BLModeFunction *constructModeFunction(
    const struct Configuration *conf) {
  struct BLModeFunction *modeFunction = 0;
  if (conf->uniform) {
    modeFunction = blModeFunctionUniformCreate(0, 1.0, 0.0);
  } else {
    modeFunction = blModeFunctionSimplifiedGaussianCreate(
        conf->waist, conf->lambda, conf->length);
  }
  return modeFunction;
}
