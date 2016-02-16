[![Build Status](https://travis-ci.org/d-meiser/BeamLaser.svg?branch=master)](https://travis-ci.org/d-meiser/BeamLaser)
[![Coverage Status](https://coveralls.io/repos/d-meiser/BeamLaser/badge.svg?branch=master&service=github)](https://coveralls.io/github/d-meiser/BeamLaser?branch=master)


# BeamLaser: Simple simulation of lasers.

The purpose of the BeamLaser code is to enable simulations
involving many independent particles interacting with few cavity
fields.  The cavity fields are treated semi-classically with
field amplitudes subject to decay and classical noise.  Atoms are
described by classical position and velocity and a quantum
mechanical state vector.  Decay and fluctuations can be
incorporated via Monte-Carlo quantum trajectory methods.  The
classical description of the atomic motion precludes simulation
of atoms in the sub-recoil regime.

The overall design of BeamLaser emphasizes pragmatism,
extensibility, and performance.


## Getting started.

Use the following commands to obtain the BeamLaser package:

```
git clone https://github.com/d-meiser/BeamLaser
cd BeamLaser
./utilities/get_sprng.sh
mkdir build
cd build
cmake ..
make -j8
make test
```

The example executables in `example` provide a good starting
point.  Better documentation is under development.


## Prerequisites

BeamLaser is being developed and tested on linux.  We use the
[centos 6 distribution](https://www.centos.org/) but any linux
distribution should work.  The following prerequisites can be
downloaded and built using the corresponding scripts in the
`utilities` directories.

- _cmake_: The [cmake](https://cmake.org/) software is needed for
  configuration.  No download script is provided for `cmake`.
  Please use your distribution's package management systems or
  download the `cmake` source from the above URL.
- _sprng_: BeamLaser minimally depends on the
  [sprng](http://www.sprng.org/Version2.0/index.html) package for
  the generation of pseudo random numbers.  We use version 2
  because it has a C interface and BeamLaser is written in C.
- _MPI_ (optional): BeamLaser can execute in parallel using MPI.
  The parallel execution model is very simplistic, so any recent
  MPI implementation should work.  Most testing has been carried
  out with [mpich](https://www.mpich.org/).
- _CGREEN_ (optional): The unit tests of BeamLaser are written
  using the [Cgreen](https://github.com/cgreen-devs/cgreen) unit
  testing framework.


### Expert configuration for hacking BeamLaser

For debugging and general development we recommend the following
configuration:

```
git clone https://github.com/d-meiser/BeamLaser
cd BeamLaser
./utilities/get_mpich.sh
./utilities/get_cgreen.sh
WITH_MPI=ON ./utilities/get_sprng.sh
mkdir build
cd build
cmake \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCGREEN_WITH_CXX=OFF \
  -DCGREEN_WITH_UNIT_TESTS=OFF \
  -DBL_BUILD_TESTS=ON \
  -DBL_ENABLE_COVERAGE=ON \
  -DBL_WITH_PEDANTIC_FLAGS=OFF \
  -DBL_WITH_WALL_FLAGS=ON \
  -DBL_WITH_MPI=ON \
  -DMPI_C_COMPILER=`pwd`/../mpich/bin/mpicc \
  -DMPI_C_INCLUDE_PATH=`pwd`/../mpich/include \
  -DMPIEXEC=`pwd`/../mpich/bin/mpirun \
  ..
make -j8
make test
```

For maximum performance change the `cmake` configuration step to:
```
cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER_FLAGS='-ffast-math -O3' \
  -DCGREEN_WITH_CXX=OFF \
  -DCGREEN_WITH_UNIT_TESTS=OFF \
  -DBL_BUILD_TESTS=ON \
  -DBL_ENABLE_COVERAGE=ON \
  -DBL_WITH_PEDANTIC_FLAGS=OFF \
  -DBL_WITH_WALL_FLAGS=ON \
  -DBL_WITH_MPI=ON \
  -DMPI_C_COMPILER=`pwd`/../mpich/bin/mpicc \
  -DMPI_C_INCLUDE_PATH=`pwd`/../mpich/include \
  -DMPIEXEC=`pwd`/../mpich/bin/mpirun \
  ..
```

## License

BeamLaser is licensed under the GPLv3.  See LICENSE for details.
