#! /usr/bin/bash

CONFIG_FILE="config.json"
PATHS_FILE="analysis_paths.json"

PROCESS=$(jq -r '.process' $CONFIG_FILE)
DISCRIMINATOR=$(jq -r '.discriminator' $CONFIG_FILE)

LOGDIR=$(jq -r --arg proc "$PROCESS" --argjson disc "$DISCRIMINATOR" '.[$proc] + "/roostats_results/out_D" + ($disc * 1000 | tostring)' $PATHS_FILE)
mkdir -p ${LOGDIR}

nohup root -q 'roostats_limits_run.cc("config.json")' > ${LOGDIR}/roostats_out.txt & 