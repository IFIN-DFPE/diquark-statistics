#! /usr/bin/bash

g++ -O3 roofit_stats.cc `root-config --cflags` -o roofit_stats -L$ROOTSYS/lib -lRooFitCore -lRooFit -lRooStats `root-config --glibs`
nohup ./roofit_stats stats.ini > results/mChi2/roofit_results/out_D925/roofit_outi.txt &