#!/bin/bash

# Usage:
# Need to link mail command to GP2 service account + Postfix
# chrom start, chrom stop, run file #, run_trial#
# ./run_pipeline.sh EUR CS 1 5 1 3

set -euo pipefail  # All processes need to succeed to keep going sequentially

# --- Run in background and disown if not already backgrounded ---
if [[ "$0" == "$BASH_SOURCE" && -z "${BACKGROUND_RUN:-}" ]]; then
  echo "[INFO] Re-launching script in background with disown..."
  BACKGROUND_RUN=1 nohup "$0" "$@" > run_pipeline.log 2>&1 & disown
  exit 0
fi

# --- Input Arguments ---
BASE_DIR=$(pwd)/google_cloud
anc=$1
analysis=$2
chrom_start=$3x
chrom_end=$4
run_script=$5
trial_run=$6

# --- Directories and Paths ---
NEXTFLOW_SCRIPT="${BASE_DIR}/long-gwas-pipeline/main.nf"
PARAMS_FILE="${BASE_DIR}/run_files/${anc}_${analysis}_chroms_tmp_dir${run_script}.yml"
RESULTS_BASE="${BASE_DIR}/run_trial${trial_run}/files/longGWAS_pipeline/results"
CHROM_DIR="CHR${chrom_start}_${chrom_end}"
DEST_BUCKET="gs://gp2_working_eu/nicole/JAN/${anc}/${analysis}/${CHROM_DIR}"
LOCAL_RESULTS_DIR="${RESULTS_BASE}/${anc}_${analysis}_CHROMS_${CHROM_DIR}"

# --- Run Nextflow Pipeline ---
echo "Running Nextflow pipeline..."
nextflow run "$NEXTFLOW_SCRIPT" -params-file "$PARAMS_FILE" -profile standard

# --- Run in background and disown if not already backgrounded ---
if [[ "$0" == "$BASH_SOURCE" && -z "$BACKGROUND_RUN" ]]; then
  echo "[INFO] Re-launching script in background with disown..."
  BACKGROUND_RUN=1 nohup "$0" "$@" > run_pipeline.log 2>&1 & disown
  exit 0
fi

# --- Transfer Results to Google Cloud Bucket ---
echo "Transferring results to GCS..."
gsutil -m cp -r "$LOCAL_RESULTS_DIR" "$DEST_BUCKET"

# --- Send Completion Email ---
EMAIL="nicole@datatecnica.com"
SUBJECT="Nextflow GWAS Pipeline Completed"
BODY="The Nextflow pipeline for ${anc}_${analysis}_CHROMS_${CHROM_DIR} has completed and results are uploaded to:\n${DEST_BUCKET}"

echo -e "$BODY" | mail -s "$SUBJECT" "$EMAIL"

echo "Done."