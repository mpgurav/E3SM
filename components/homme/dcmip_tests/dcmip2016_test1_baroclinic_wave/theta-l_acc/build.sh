#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make theta-l_acc-nlev30"
  make -j theta-l_acc-nlev30
cd $cwd
