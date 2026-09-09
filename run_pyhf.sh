#!/bin/bash
set -euo pipefail

MASTER_LOG="master_run.txt"
: > "$MASTER_LOG"

# keep a master log in the main directory
exec > >(tee -a "$MASTER_LOG") 2>&1

echo "=========================================================="
echo "Starting master run at $(date)"
echo "=========================================================="

# list the configs you want to run
# use your actual config filenames here
CONFIGS=(
    # "configs/ChiChi/cls_as_yuu02_mChi20_D900.json"
    # "configs/ChiChi/cls_tb_yuu02_mChi20_D900.json"    
    
    "configs/ChiChi/ul_as_yuu04_mChi20_D900.json"
    # "configs/ChiChi/ul_tb_yuu04_mChi20_D900.json"
    "configs/ChiChi/ul_as_yuu02_mChi20_D900.json"
    # "configs/ChiChi/ul_tb_yuu02_mChi20_D900.json"

    # "configs/ChiChi/pval_as_yuu02_mChi20_D900.json"
    # "configs/ChiChi/pval_tb_yuu02_mChi20_D900.json"
)

for cfg in "${CONFIGS[@]}"; do
    if [ ! -f "$cfg" ]; then
        echo "Skipping missing config: $cfg"
        continue
    fi

    # extract process from the config
    PROCESS=$(jq -r '.process' "$cfg")
    RESULT_DIR=$(jq -r --arg proc "$PROCESS" '.[$proc]' analysis_paths.json)/pyhf_output

    mkdir -p "$RESULT_DIR"

    # log file stored next to the results, using only the config filename
    CFG_NAME=$(basename "$cfg")
    LOG_FILE="${RESULT_DIR}/${CFG_NAME%.json}.log"

    echo ""
    echo "----------------------------------------------------------"
    echo "Running: $cfg at $(date)"
    echo "Result directory: $RESULT_DIR"
    echo "Log file: $LOG_FILE"
    echo "----------------------------------------------------------"

    python3 run_stats.py --config "$cfg" --paths analysis_paths.json > "$LOG_FILE" 2>&1

    echo "Finished $cfg at $(date)"
done

echo ""
echo "=========================================================="
echo "All runs completed at $(date)"
echo "=========================================================="