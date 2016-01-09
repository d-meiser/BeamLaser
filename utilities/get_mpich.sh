#!/bin/sh

if [ -f mpich/include/mpi.h ]; then
  echo "mpich/include/mpi.h found."
fi
if [ -f mpich/lib/libmpich.so ]; then
  echo "libmpich.so found -- nothing to build."
else
  if [ -f mpich-3.2/src/mpi/init/mpi_init.h ]; then
    echo "mpich sources found -- don't need to download"
  else
    echo "Downloading mpich source."
    wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
    tar xfz mpich-3.2.tar.gz
    rm mpich-3.2.tar.gz
  fi
  echo "configuring and building mpich."
  cd mpich-3.2
  ./configure \
          --prefix=`pwd`/../mpich \
          --enable-static=false \
          --enable-alloca=true \
          --disable-long-double \
          --enable-threads=single \
          --enable-fortran=no \
          --enable-fast=all \
          --enable-g=none \
          --enable-timing=none
  make -j4
  make install
  cd -
fi

