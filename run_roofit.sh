#! /usr/bin/bash

set -e
g++ -O3 -fopenmp roofit_stats.cc `root-config --cflags` -o roofit_stats -L$ROOTSYS/lib -lRooFitCore -lRooFit -lRooStats `root-config --glibs`
export OMP_NUM_THREADS=8
./roofit_stats stats.ini
# nohup ./roofit_stats stats.ini > results/mChi2/roofit_results/out_D925/roofit_out_S825.txt &