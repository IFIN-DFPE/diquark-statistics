#!/usr/bin/bash

# Exit immediately if a command crashes
set -e

# ==============================================================================
# 0. SETUP GLOBAL LOGGING
# This line captures ALL standard output and errors from this bash script
# and saves it to 'master_run_log.txt', while also printing it to your screen.
# ==============================================================================
exec > >(tee -i master_run_log.txt) 2>&1

echo "=========================================================="
echo "Starting analysis runs at $(date)"
echo "=========================================================="

# ==============================================================================
# 1. COMPILE ALL PROGRAMS UPFRONT
# ==============================================================================
ROOT_FLAGS=$(root-config --cflags)
ROOT_LIBS="-L$ROOTSYS/lib -lRooFitCore -lRooFit -lRooStats $(root-config --glibs)"

echo "Compiling roofit_pval..."
g++ -O3 -g roofit_pval.cc $ROOT_FLAGS -o roofit_pval $ROOT_LIBS

echo "Compiling roofit_stats..."
g++ -O3 -g roofit_stats.cc $ROOT_FLAGS -o roofit_stats $ROOT_LIBS

echo "Compiling roofit_pval_pdf_halved..."
g++ -O3 -g roofit_pval_pdf_halved.cc $ROOT_FLAGS -o roofit_pval_pdf_halved $ROOT_LIBS

echo "Compiling roofit_pval_pdf_quarter..."
g++ -O3 -g roofit_pval_pdf_quarter.cc $ROOT_FLAGS -o roofit_pval_pdf_quarter $ROOT_LIBS

echo "All compilations successful."
echo ""


# ==============================================================================
# 2A. FUNCTION FOR COMPILED EXECUTABLES (roofit_results)
# ==============================================================================
PATHS_FILE="analysis_paths.json"

run_analysis() {
    local EXE_NAME=$1
    local CONFIG_FILE=$2

    echo "----------------------------------------------------------"
    echo "Running: ./$EXE_NAME with $CONFIG_FILE at $(date)" 
    echo "----------------------------------------------------------"

    # Verify config file exists
    if [ ! -f "$CONFIG_FILE" ]; then
        echo "Warning: $CONFIG_FILE not found. Skipping..."
        return 0 # Skips to the next command without crashing the whole script
    fi

    # Extract parameters
    PROCESS=$(jq -r '.process' "$CONFIG_FILE")
    DISCRIMINATOR=$(jq -r '.discriminator' "$CONFIG_FILE")

    # Create log directory
    LOGDIR=$(jq -r --arg proc "$PROCESS" --argjson disc "$DISCRIMINATOR" '.[$proc] + "/roofit_results/out_D" + ($disc * 1000 | tostring)' "$PATHS_FILE")
    mkdir -p "${LOGDIR}"

    # Name the log file dynamically based on the executable name 
    # (so roofit_stats doesn't overwrite roofit_pval's log)
    LOGFILE="${LOGDIR}/${EXE_NAME}_out.txt"
    echo "Program logs will be saved to: $LOGFILE"

    # Execute synchronously
    ./$EXE_NAME "$CONFIG_FILE" > "$LOGFILE" 2>&1

    echo "Finished ./$EXE_NAME with $CONFIG_FILE at $(date)."
    echo ""
}

# ==============================================================================
# 2B. FUNCTION FOR ROOT MACROS (roostats_results)
# ==============================================================================
run_roostats_macro() {
    local MACRO_NAME=$1
    local CONFIG_FILE=$2

    echo "----------------------------------------------------------"
    echo "Running ROOT Macro: $MACRO_NAME with $CONFIG_FILE at $(date)"
    echo "----------------------------------------------------------"

    if [ ! -f "$CONFIG_FILE" ]; then
        echo "Warning: $CONFIG_FILE not found. Skipping..."
        return 0
    fi

    PROCESS=$(jq -r '.process' "$CONFIG_FILE")
    DISCRIMINATOR=$(jq -r '.discriminator' "$CONFIG_FILE")

    # Creates roostats_results directory (Notice the change here!)
    LOGDIR=$(jq -r --arg proc "$PROCESS" --argjson disc "$DISCRIMINATOR" '.[$proc] + "/roostats_results/out_D" + ($disc * 1000 | tostring)' "$PATHS_FILE")
    mkdir -p "${LOGDIR}"

    # Strip the .cc extension to make a clean log file name
    local MACRO_BASE="${MACRO_NAME%.cc}"
    LOGFILE="${LOGDIR}/${MACRO_BASE}_out.txt"
    echo "Macro logs will be saved to: $LOGFILE"

    # Execute ROOT macro synchronously in batch mode (-b)
    # The escaped quotes (\" \") ensure the config file name is passed as a string to C++
    root -q -b "${MACRO_NAME}(\"${CONFIG_FILE}\")" > "$LOGFILE" 2>&1

    echo "Finished ROOT Macro $MACRO_NAME with $CONFIG_FILE at $(date)."
    echo ""
}



# ==============================================================================
# 3. SCHEDULE YOUR RUNS HERE
# Call the functions: 
#           run_analysis "program_name" "config_file.json"
#           run_roostats_macro "macro_name.cc" "config_file.json"
# They will run strictly sequentially, waiting for the previous to finish.
# ==============================================================================

# run_roostats_macro "roostats_limits_run.cc" "ChiChi/configs/config_pval_yuu02_mChi15_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "ChiChi/configs/config_pval_yuu02_mChi20_D800.json"
# run_roostats_macro "roostats_limits_run.cc" "ChiChi/configs/config_pval_yuu02_mChi20_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "ChiChi/configs/config_pval_yuu02_mChi20_D925.json"
# run_roostats_macro "roostats_limits_run.cc" "ChiChi/configs/config_pval_yuu04_mChi15_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "ChiChi/configs/config_pval_yuu04_mChi20_D900.json"

# run_analysis "roofit_pval" "ChiChi/configs/config_pval_yuu02_mChi15_D900.json"
# run_analysis "roofit_pval" "ChiChi/configs/config_pval_yuu02_mChi20_D800.json"
# run_analysis "roofit_pval" "ChiChi/configs/config_pval_yuu02_mChi20_D900.json"
# run_analysis "roofit_pval" "ChiChi/configs/config_pval_yuu02_mChi20_D900_v1.json"
# run_analysis "roofit_pval" "ChiChi/configs/config_pval_yuu02_mChi20_D900_v2.json"
# run_analysis "roofit_pval" "ChiChi/configs/config_pval_yuu02_mChi20_D925.json"
# run_analysis "roofit_pval" "ChiChi/configs/config_pval_yuu04_mChi20_D900.json"

# run_analysis "roofit_pval_pdf_halved" "ChiChi/configs/config_pval_yuu02_mChi20_D900.json"
# run_analysis "roofit_pval_pdf_halved" "ChiChi/configs/config_pval_yuu04_mChi20_D900.json"

# run_analysis "roofit_stats" "ChiChi/configs/config_cls_yuu02_mChi15_D900.json"
# run_analysis "roofit_stats" "ChiChi/configs/config_cls_yuu02_mChi20_D900.json"





# run_roostats_macro "roostats_limits_run.cc" "uChi/configs/config_pval_yuu02_mChi15_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "uChi/configs/config_pval_yuu02_mChi20_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "uChi/configs/config_pval_yuu04_mChi15_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "uChi/configs/config_pval_yuu04_mChi20_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "uChi/configs/config_pval_yuChi05_mChi20_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "uChi/configs/config_pval_uZt_mChi20_D900.json"
# run_roostats_macro "roostats_limits_run.cc" "uChi/configs/config_pval_ubbbart_mChi20_D900.json"

# run_analysis "roofit_pval" "uChi/configs/config_pval_yuu02_mChi15_D900.json"
# run_analysis "roofit_pval" "uChi/configs/config_pval_yuu02_mChi20_D900.json"
# run_analysis "roofit_pval" "uChi/configs/config_pval_yuu04_mChi20_D900.json"
# run_analysis "roofit_pval" "uChi/configs/config_pval_yuChi05_mChi20_D900.json"
# run_analysis "roofit_pval" "uChi/configs/config_pval_uZt_mChi20_D900.json"
# run_analysis "roofit_pval" "uChi/configs/config_pval_ubbbart_mChi20_D900.json"

# run_analysis "roofit_pval_pdf_halved" "uChi/configs/config_pval_yuu02_mChi20_D900.json"
# run_analysis "roofit_pval_pdf_halved" "uChi/configs/config_pval_yuu04_mChi20_D900.json"
# run_analysis "roofit_pval_pdf_halved" "uChi/configs/config_pval_yuChi05_mChi20_D900.json"

# run_analysis "roofit_pval_pdf_quarter" "uChi/configs/config_pval_yuu02_mChi20_D900.json"
# run_analysis "roofit_pval_pdf_quarter" "uChi/configs/config_pval_yuu04_mChi20_D900.json"
# run_analysis "roofit_pval_pdf_quarter" "uChi/configs/config_pval_yuChi05_mChi20_D900.json"

run_analysis "roofit_stats" "uChi/configs/config_cls_yuu02_mChi15_D900.json"
run_analysis "roofit_stats" "uChi/configs/config_cls_yuu02_mChi20_D900.json"
run_analysis "roofit_stats" "uChi/configs/config_cls_ubbbart_mChi20_D900.json"
run_analysis "roofit_stats" "uChi/configs/config_cls_uZt_mChi20_D900.json"



echo "=========================================================="
echo "All sequential runs completed successfully at $(date)!"
echo "=========================================================="