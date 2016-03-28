#!/bin/sh

if [ -f guile/bin/guile ]; then
    echo "guile found -- nothing to build."
else
  if [ -f guile-2.0.11/libguile.h ]; then
    echo "guile sources found -- don't need to download"
  else
    echo "Downloading guile source."
    wget ftp://ftp.gnu.org/gnu/guile/guile-2.0.11.tar.gz
    tar xf guile-2.0.11.tar.gz
    rm guile-2.0.11.tar.gz
  fi
  echo "configuring and building guile."
  cd guile-2.0.11
  ./configure \
          --prefix=`pwd`/../guile
  make -j4
  make install
  cd -
fi

