#!/bin/sh

if [ ! -e sprng5 ]; then
  echo "Downloading sprng source."
  wget http://www.sprng.org/Version5.0/sprng5.tar.bz2
  tar xjf sprng5.tar.bz2
  rm sprng5.tar.bz2
fi

cd sprng5
./configure
make
cd -
