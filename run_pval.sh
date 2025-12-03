#! /usr/bin/bash

set -e

g++ -O3 -g roofit_pval.cc `root-config --cflags` -o roofit_pval -L$ROOTSYS/lib -lRooFitCore -lRooFit -lRooStats `root-config --glibs`

CONFIG_FILE="config.json"
PATHS_FILE="analysis_paths.json"

PROCESS=$(jq -r '.process' $CONFIG_FILE)
DISCRIMINATOR=$(jq -r '.discriminator' $CONFIG_FILE)

LOGDIR=$(jq -r --arg proc "$PROCESS" --argjson disc "$DISCRIMINATOR" '.[$proc] + "/roofit_results/out_D" + ($disc * 1000 | tostring)' $PATHS_FILE)
mkdir -p ${LOGDIR}

nohup ./roofit_pval config.json > ${LOGDIR}/roofit_pval_out.txt &
