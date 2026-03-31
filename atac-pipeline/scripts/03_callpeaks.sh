#!/bin/bash
set -euo pipefail

#######################################################################################
##                         ATAC-seq PEAK CALLING WITH MACS3                          ##
#######################################################################################


#### OVERVIEW ####
# Calls peaks on deduplicated BAM files using MACS3.
# Expects clean BAMs produced by 02_atac_preprocess_align.sh.


#### INPUT DATA ####
# Samples from mESCs with Ascl1-GR-HA
# 3 conditions: mESC, EpiLC, NE
# 3 timepoints per condition: 0h, 6h, 24h after dexamethasone induction of Ascl1
# Sample sheet: first column is the full sample name e.g. C11_EpiLC_HA_1
# Clean BAMs expected at: <MAINDIR>/cleanbams/<sample>.nodup.bam


#### USAGE ####
# bash 03_call_peaks.sh <main_directory> <sample_sheet> [parallel_jobs]
#
# Arguments:
#   $1  MAINDIR       - path to the main project directory
#   $2  SAMPLESHEET   - path to CSV sample sheet (first column = sample name, with header)
#   $3  PARALLEL_JOBS - number of samples to process in parallel (default: 4)
#
# Example:
#   bash 03_call_peaks.sh /path/to/project /path/to/samplesheet.csv 8


#### ARGUMENT HANDLING ####
if [[ $# -lt 2 ]]; then
    echo "Usage: bash 03_call_peaks.sh <main_directory> <sample_sheet> [parallel_jobs]"
    exit 1
fi

MAINDIR=$1
SAMPLESHEET=$2
PARALLEL_JOBS=${3:-4}

# Genome size for mm39 (Mus musculus)
GENOME_SIZE=2468088461

FILTERDIR=$MAINDIR/cleanbams
PEAKDIR=$MAINDIR/macs3

echo "Main directory : $MAINDIR"
echo "Sample sheet   : $SAMPLESHEET"
echo "Parallel jobs  : $PARALLEL_JOBS"
echo "Genome size    : $GENOME_SIZE (mm39)"
echo "Clean BAMs     : $FILTERDIR"
echo "Output dir     : $PEAKDIR"


#### CHECK CLEAN BAMS DIRECTORY EXISTS ####
if [[ ! -d "$FILTERDIR" ]]; then
    echo "[ERROR]: Clean BAMs directory not found at $FILTERDIR."
    echo "         Please run 02_atac_preprocess_align.sh first."
    exit 1
fi

mkdir -p "$PEAKDIR"


#### GET SAMPLE NAMES ####
SAMPLES=()
while IFS=$'\r\n' read -r line
do
    VALUES=(${line//,/ })
    SAMPLES+=("${VALUES[0]}")
    echo "Sample: ${VALUES[0]}"
done <<< $(tail +2 "$SAMPLESHEET")


#######################################################################################
##                           Call peaks with MACS3                                   ##
#######################################################################################

echo "Calling peaks with MACS3..."

callPeaks_macs3() {
    local s="$1"
    local input_bam="$FILTERDIR/${s}.nodup.bam"
    local output_peak="$PEAKDIR/${s}_peaks.narrowPeak"

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: $input_bam not found. Skipping sample $s."
        return
    fi

    if [[ -f "$output_peak" ]]; then
        echo "Peaks for sample $s already called. Skipping..."
        return
    fi

    echo "Calling peaks for sample $s..."
    macs3 callpeak \
        -t "$input_bam" \
        -f BAMPE \
        -n "$s" \
        -g "$GENOME_SIZE" \
        --outdir "$PEAKDIR" \
        -B
    echo "Peak calling complete for sample $s"
}

export -f callPeaks_macs3
export FILTERDIR PEAKDIR GENOME_SIZE

printf "%s\n" "${SAMPLES[@]}" | xargs -I {} -P "$PARALLEL_JOBS" bash -c 'callPeaks_macs3 "$@"' _ {}
wait

echo "Peak calling complete. Output in: $PEAKDIR"
