#!/bin/bash
set -euo pipefail

#######################################################################################
##                    ATAC-seq BIGWIG GENERATION WITH TMM SCALING                    ##
#######################################################################################


#### OVERVIEW ####
# Generates TMM-normalised bigwig files from deduplicated BAM files using bamCoverage.
# Scale factors are read from a CSV produced by the DiffBind TMM normalisation step
# in R (see companion R script). Each sample's bamCoverage scale factor (1/TMM factor)
# is applied so that bigwigs are comparable across samples.


#### INPUT DATA ####
# Samples from mESCs with Ascl1-GR-HA
# 3 conditions: mESC, EpiLC, NE
# 3 timepoints per condition: 0h, 6h, 24h after dexamethasone induction of Ascl1
# BAM files expected as:   <BAMDIR>/<sample>.nodup.bam
# Scale factor CSV expected columns: sample, bam_coverage_scale
# (produced by DiffBind TMM normalisation — see R script)


#### USAGE ####
# bash 04_bigwig_tmm.sh <bam_directory> <scale_factor_csv> [threads] [bin_size]
#
# Arguments:
#   $1  BAMDIR           - directory containing deduplicated .nodup.bam files
#   $2  SCALE_FACTOR_CSV - path to CSV with columns: sample, bam_coverage_scale
#   $3  THREADS          - number of threads per bamCoverage job (default: 8)
#   $4  BINSIZE          - bigwig bin size in bp (default: 10)
#
# Example:
#   bash 04_bigwig_tmm.sh /path/to/bam /path/to/TMM_scale_factors.csv 8 10


#### ARGUMENT HANDLING ####
if [[ $# -lt 2 ]]; then
    echo "Usage: bash 04_bigwig_tmm.sh <bam_directory> <scale_factor_csv> [threads] [bin_size]"
    exit 1
fi

BAMDIR=$1
SCALE_FACTOR_CSV=$2
THREADS=${3:-8}
BINSIZE=${4:-10}

# Output bigwig directory sits alongside the bam directory
BIGWIGDIR=$(dirname "$BAMDIR")/bigwig
mkdir -p "$BIGWIGDIR"

echo "BAM directory      : $BAMDIR"
echo "Scale factor CSV   : $SCALE_FACTOR_CSV"
echo "Bigwig directory   : $BIGWIGDIR"
echo "Threads            : $THREADS"
echo "Bin size           : $BINSIZE bp"


#### CHECK INPUTS ####
if [[ ! -d "$BAMDIR" ]]; then
    echo "[ERROR]: BAM directory not found at $BAMDIR."
    exit 1
fi

if [[ ! -f "$SCALE_FACTOR_CSV" ]]; then
    echo "[ERROR]: Scale factor CSV not found at $SCALE_FACTOR_CSV."
    exit 1
fi


#######################################################################################
##                        Generate TMM-scaled bigwigs                                ##
#######################################################################################

echo "Generating TMM-scaled bigwigs..."

# Skip header line and process each sample in the CSV
tail -n +2 "$SCALE_FACTOR_CSV" | while IFS=',' read -r sample lib_size tmm_scale_factor bam_coverage_scale
do
    # Strip any quotes or whitespace from CSV fields
    sample=$(echo "$sample" | tr -d '"' | xargs)
    bam_coverage_scale=$(echo "$bam_coverage_scale" | tr -d '"' | xargs)

    input_bam="$BAMDIR/${sample}.nodup.bam"
    output_bw="$BIGWIGDIR/${sample}.TMM.bw"

    # Check BAM exists
    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: BAM file not found for sample $sample at $input_bam. Skipping..."
        continue
    fi

    # Check BAI index exists
    if [[ ! -f "${input_bam}.bai" && ! -f "${input_bam%.bam}.bai" ]]; then
        echo "[ERROR]: BAM index not found for $input_bam. Skipping..."
        continue
    fi

    # Skip if bigwig already exists
    if [[ -f "$output_bw" ]]; then
        echo "Bigwig for sample $sample already exists. Skipping..."
        continue
    fi

    echo "Processing $sample (scale factor: $bam_coverage_scale)..."
    bamCoverage \
        --bam "$input_bam" \
        --outFileName "$output_bw" \
        --outFileFormat bigwig \
        --scaleFactor "$bam_coverage_scale" \
        --normalizeUsing None \
        --binSize "$BINSIZE" \
        --numberOfProcessors "$THREADS" \
        --extendReads

    echo "Bigwig complete for sample $sample"
done

echo "All bigwigs written to: $BIGWIGDIR"
