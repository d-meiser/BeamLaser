#!/bin/sh

if [ -f sprng2.0/lib/libsprng.a ]; then
  echo "libsprng.a found -- nothing  to build."
else
  if [ -f sprng2.0/SRC/sprng.h ]; then
    echo "Already have sprng source -- no need to download."
  else
    echo "Downloading sprng source."
    wget http://www.sprng.org/Version2.0/sprng2.0b.tar.gz
    tar xfz sprng2.0b.tar.gz
    rm sprng2.0b.tar.gz
  fi
  echo "building libsprng.a"
  cd sprng2.0
  sed -i.bak s/^PMLCGDEF/#PMLCGDEF/ make.CHOICES
  sed -i.bak s/^GMPLIB/#GMPLIB/ make.CHOICES
  sed -i.bak s/g77/echo/ SRC/make.INTEL
  make src
  cd -
fi

