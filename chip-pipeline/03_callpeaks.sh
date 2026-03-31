#!/bin/bash
set -euo pipefail

#######################################################################################
##                         ChIP-seq PEAK CALLING WITH MACS3                          ##
#######################################################################################


#### OVERVIEW ####
# Calls peaks on duplicate-marked BAM files using MACS3 with matched input controls.
# Expects clean BAMs produced by chip_02_preprocess_align.sh.
# A CSV input-matching file specifies which input control to use for each ChIP sample.


#### INPUT DATA ####
# Samples from HA-Ascl1 ChIP in three conditions: mESC, EpiLC, NE
# Clean BAMs expected at: <MAINDIR>/cleanbams/<sample>_markdup.bam
# Input matching CSV expected columns: sample, target, input
#   sample - output peak name
#   target - ChIP BAM sample name
#   input  - matched input control BAM sample name


#### USAGE ####
# bash chip_03_call_peaks.sh <main_directory> <input_matching_csv> [parallel_jobs]
#
# Arguments:
#   $1  MAINDIR           - path to the main project directory
#   $2  INPUT_MATCHING    - path to CSV with columns: sample, target, input
#   $3  PARALLEL_JOBS     - number of samples to process in parallel (default: 4)
#
# Example:
#   bash chip_03_call_peaks.sh /path/to/project /path/to/chip_inputmatching.csv 8


#### ARGUMENT HANDLING ####
if [[ $# -lt 2 ]]; then
    echo "Usage: bash chip_03_call_peaks.sh <main_directory> <input_matching_csv> [parallel_jobs]"
    exit 1
fi

MAINDIR=$1
INPUT_MATCHING=$2
PARALLEL_JOBS=${3:-4}

# Genome size for mm39 (Mus musculus)
# Source: https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
GENOME_SIZE=2468088461

FILTERDIR=$MAINDIR/cleanbams
PEAKDIR=$MAINDIR/macs3

echo "Main directory  : $MAINDIR"
echo "Input matching  : $INPUT_MATCHING"
echo "Parallel jobs   : $PARALLEL_JOBS"
echo "Genome size     : $GENOME_SIZE (mm39)"
echo "Clean BAMs      : $FILTERDIR"
echo "Output dir      : $PEAKDIR"


#### CHECK INPUTS ####
if [[ ! -d "$FILTERDIR" ]]; then
    echo "[ERROR]: Clean BAMs directory not found at $FILTERDIR."
    echo "         Please run chip_02_preprocess_align.sh first."
    exit 1
fi

if [[ ! -f "$INPUT_MATCHING" ]]; then
    echo "[ERROR]: Input matching CSV not found at $INPUT_MATCHING."
    exit 1
fi

mkdir -p "$PEAKDIR"


#######################################################################################
##                           Call peaks with MACS3                                   ##
#######################################################################################

echo "Calling peaks with MACS3..."

run_macs3() {
    local sample="$1"
    local target="$2"
    local input="$3"

    local target_bam="$FILTERDIR/${target}_markdup.bam"
    local input_bam="$FILTERDIR/${input}_markdup.bam"
    local output_peak="$PEAKDIR/${sample}_peaks.narrowPeak"

    if [[ ! -f "$target_bam" ]]; then
        echo "[ERROR]: ChIP BAM not found for $target at $target_bam. Skipping $sample."
        return
    fi

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: Input BAM not found for $input at $input_bam. Skipping $sample."
        return
    fi

    if [[ -f "$output_peak" ]]; then
        echo "Peaks for $sample already called. Skipping..."
        return
    fi

    echo "Calling peaks for $sample (ChIP: $target, input: $input)..."
    macs3 callpeak \
        -t "$target_bam" \
        -c "$input_bam" \
        -g "$GENOME_SIZE" \
        -f BAMPE \
        -n "$sample" \
        --outdir "$PEAKDIR"

    echo "Peak calling complete for $sample"
}

export -f run_macs3
export FILTERDIR PEAKDIR GENOME_SIZE

# Read input matching CSV and dispatch jobs in parallel
# CSV format: sample,target,input (with header)
tail -n +2 "$INPUT_MATCHING" | while IFS=',' read -r sample target input rest; do
    sample=$(echo "$sample" | tr -d '"' | xargs)
    target=$(echo "$target" | tr -d '"' | xargs)
    input=$(echo "$input"  | tr -d '"' | xargs)
    echo "$sample $target $input"
done | xargs -P "$PARALLEL_JOBS" -n 3 bash -c 'run_macs3 "$@"' _

echo "Peak calling complete. Output in: $PEAKDIR"
