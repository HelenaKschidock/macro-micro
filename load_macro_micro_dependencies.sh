#!/usr/bin/env bash

PRECICE_PREFIX=~/software/precice # set this to your selected prefix
export PATH=$PRECICE_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PRECICE_PREFIX/lib:$LD_LIBRARY_PATH
export CPATH=$PRECICE_PREFIX/include:$CPATH
export PKG_CONFIG_PATH=$PRECICE_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
export CMAKE_PREFIX_PATH=$PRECICE_PREFIX:$CMAKE_PREFIX_PATH
module use /usr/local.nfs/sgs/modulefiles
module load ub2004/boost/1.75.0
module load ub2004/libxml2/2.9.10
module load cmake/3.18.2
