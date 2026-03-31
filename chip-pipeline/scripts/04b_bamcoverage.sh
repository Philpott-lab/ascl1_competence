#!/bin/bash
set -euo pipefail

#######################################################################################
##               ChIP-seq BIGWIG GENERATION AND REPLICATE MERGING                   ##
#######################################################################################


#### OVERVIEW ####
# Generates RLE-normalised bigwig files from duplicate-marked BAM files using
# bamCoverage, then merges replicates per condition using bigwigAverage.
# Scale factors are read from the text file produced by chip_05a_diffbind_count_normalise.R.


#### INPUT DATA ####
# Samples from HA-Ascl1 ChIP in three conditions: mESC, EpiLC, NE
# BAM files expected as:   <BAMDIR>/<sample>_markdup.bam
# Scale factor file expected columns (space-delimited): bamID NormFacs InvNormFacs
# (produced by chip_05a_diffbind_count_normalise.R)


#### USAGE ####
# bash chip_05b_bigwig.sh <bam_directory> <scale_factor_file> [threads] [bin_size]
#
# Arguments:
#   $1  BAMDIR           - directory containing _markdup.bam files
#   $2  SCALE_FACTOR_FILE - path to space-delimited scale factor file
#   $3  THREADS          - number of threads per bamCoverage job (default: 8)
#   $4  BINSIZE          - bigwig bin size in bp (default: 1)
#
# Example:
#   bash chip_05b_bigwig.sh /path/to/bam /path/to/ChIPseq_ASCL1_scalefactors.txt 8 1


#### ARGUMENT HANDLING ####
if [[ $# -lt 2 ]]; then
    echo "Usage: bash chip_05b_bigwig.sh <bam_directory> <scale_factor_file> [threads] [bin_size]"
    exit 1
fi

BAMDIR=$1
SCALE_FACTOR_FILE=$2
THREADS=${3:-8}
BINSIZE=${4:-1}

# Output bigwig directory sits alongside the bam directory
BIGWIGDIR=$(dirname "$BAMDIR")/bigwig
mkdir -p "$BIGWIGDIR"

echo "BAM directory      : $BAMDIR"
echo "Scale factor file  : $SCALE_FACTOR_FILE"
echo "Bigwig directory   : $BIGWIGDIR"
echo "Threads            : $THREADS"
echo "Bin size           : $BINSIZE bp"


#### CHECK INPUTS ####
if [[ ! -d "$BAMDIR" ]]; then
    echo "[ERROR]: BAM directory not found at $BAMDIR."
    exit 1
fi

if [[ ! -f "$SCALE_FACTOR_FILE" ]]; then
    echo "[ERROR]: Scale factor file not found at $SCALE_FACTOR_FILE."
    exit 1
fi


#######################################################################################
##                    1. Generate RLE-scaled bigwigs per sample                      ##
#######################################################################################

echo "Generating RLE-scaled bigwigs..."

# Skip header line; columns are: bamID NormFacs InvNormFacs
tail -n +2 "$SCALE_FACTOR_FILE" | while read -r bamID NormFacs InvNormFacs
do
    input_bam="$BAMDIR/${bamID}"
    sample="${bamID%_markdup.bam}"
    output_bw="$BIGWIGDIR/${sample}_RLE.bigwig"

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: BAM file not found at $input_bam. Skipping..."
        continue
    fi

    if [[ ! -f "${input_bam%.bam}.bai" && ! -f "${input_bam}.bai" ]]; then
        echo "[ERROR]: BAM index not found for $input_bam. Skipping..."
        continue
    fi

    if [[ -f "$output_bw" ]]; then
        echo "Bigwig for $sample already exists. Skipping..."
        continue
    fi

    echo "Processing $sample (RLE scale factor: $InvNormFacs)..."
    bamCoverage \
        --bam "$input_bam" \
        --outFileName "$output_bw" \
        --outFileFormat bigwig \
        --scaleFactor "$InvNormFacs" \
        --normalizeUsing None \
        --binSize "$BINSIZE" \
        --extendReads \
        --numberOfProcessors "$THREADS"

    echo "Bigwig complete for $sample"
done

echo "Per-sample bigwigs complete"


#######################################################################################
##                    2. Merge replicates per condition with bigwigAverage           ##
#######################################################################################

echo "Merging replicates by condition..."

# mESC replicates
if [[ -f "$BIGWIGDIR/gG5_mESC_HA_RLE_merged.bigwig" ]]; then
    echo "mESC merged bigwig already exists. Skipping..."
else
    echo "Merging mESC replicates..."
    bigwigAverage \
        -b "$BIGWIGDIR"/*mESC*RLE.bigwig \
        -o "$BIGWIGDIR/gG5_mESC_HA_RLE_merged.bigwig" \
        -bs "$BINSIZE" \
        -p "$THREADS" \
        -v
fi

# EpiLC replicates
if [[ -f "$BIGWIGDIR/gG5_EpiLC_HA_RLE_merged.bigwig" ]]; then
    echo "EpiLC merged bigwig already exists. Skipping..."
else
    echo "Merging EpiLC replicates..."
    bigwigAverage \
        -b "$BIGWIGDIR"/*EpiLC*RLE.bigwig \
        -o "$BIGWIGDIR/gG5_EpiLC_HA_RLE_merged.bigwig" \
        -bs "$BINSIZE" \
        -p "$THREADS" \
        -v
fi

# NE replicates
if [[ -f "$BIGWIGDIR/gG5_NE_HA_RLE_merged.bigwig" ]]; then
    echo "NE merged bigwig already exists. Skipping..."
else
    echo "Merging NE replicates..."
    bigwigAverage \
        -b "$BIGWIGDIR"/*NE*RLE.bigwig \
        -o "$BIGWIGDIR/gG5_NE_HA_RLE_merged.bigwig" \
        -bs "$BINSIZE" \
        -p "$THREADS" \
        -v
fi

echo "Replicate merging complete"
echo "All bigwigs written to: $BIGWIGDIR"
