#! /usr/bin/bash

nohup root -q 'roostats_limits_run.cc("results/mChi2/signal_yields/sig_bkg_D900.csv")' > results/mChi2/roostats_results/out_D900/roostats_out.txt &