#!/bin/bash
set -euo pipefail

#######################################################################################
##                    ChIP-seq DATA PRE-PROCESSING & ALIGNMENT                       ##
#######################################################################################


#### OVERVIEW ####
# 1. FastQC checks the quality of the raw fastq sequencing reads
# 2. Trim Galore trims adapters and performs post-trim QC
# 3. Reads are aligned to the combined mm39 + dm6 genome using bowtie2
# 4. mm39 and dm6 reads are separated by chromosome name prefix
# 5. Read groups are added and duplicates marked with Picard
# 6. Insert size metrics are collected


#### INPUT DATA ####
# Samples from HA-Ascl1 ChIP in three conditions: mESC, EpiLC, NE
# Sample sheet: first column is the full sample name e.g. C11_EpiLC_HA_1
# FASTQ files expected as: <MAINDIR>/fastq/<sample>_R1.fastq.gz and _R2.fastq.gz
# Combined mm39/dm6 bowtie2 index expected at: <MAINDIR>/indexes/mm39_dm6
# (run chip_01_prepare_index.sh first)


#### CONDA ENVIRONMENT ####
# conda create -n chipseq
# conda install -n chipseq \
#   cutadapt fastqc multiqc trim-galore bowtie2 samtools sambamba \
#   picard bedtools bamtools preseq pysam ucsc-bedgraphtobigwig \
#   deeptools homer subread
# See environment.yml for pinned versions.


#### USAGE ####
# bash chip_02_preprocess_align.sh <main_directory> <sample_sheet> [threads]
#
# Arguments:
#   $1  MAINDIR     - path to the main project directory
#   $2  SAMPLESHEET - path to CSV sample sheet (first column = sample name, with header)
#   $3  THREADS     - number of threads for alignment and samtools (default: 8)
#
# Example:
#   bash chip_02_preprocess_align.sh /path/to/project /path/to/samplesheet.csv 40


#### ARGUMENT HANDLING ####
if [[ $# -lt 2 ]]; then
    echo "Usage: bash chip_02_preprocess_align.sh <main_directory> <sample_sheet> [threads]"
    exit 1
fi

MAINDIR=$1
SAMPLESHEET=$2
THREADS=${3:-8}

INDEXDIR=$MAINDIR/indexes
INDEX_PREFIX="$INDEXDIR/mm39_dm6"

echo "Main directory : $MAINDIR"
echo "Sample sheet   : $SAMPLESHEET"
echo "Threads        : $THREADS"


#### CHECK BOWTIE2 INDEX EXISTS ####
if [[ ! -f "${INDEX_PREFIX}.1.bt2" ]]; then
    echo "[ERROR]: Bowtie2 index not found at ${INDEX_PREFIX}."
    echo "         Please run chip_01_prepare_index.sh first."
    exit 1
fi


#### GET SAMPLE NAMES ####
SAMPLES=()
while IFS=$'\r\n' read -r line
do
    VALUES=(${line//,/ })
    SAMPLES+=("${VALUES[0]}")
    echo "Sample: ${VALUES[0]}"
done <<< $(tail +2 "$SAMPLESHEET")


#######################################################################################
##                          1. Initial QC with FastQC                                ##
#######################################################################################

READDIR=$MAINDIR/fastq
QCDIR=$MAINDIR/fastqc
mkdir -p "$QCDIR"

echo "Performing initial FastQC..."

run_fastqc() {
    local s="$1"
    local r1_file="$READDIR/${s}_R1.fastq.gz"
    local r2_file="$READDIR/${s}_R2.fastq.gz"
    local output_r1="$QCDIR/${s}_R1_fastqc.zip"
    local output_r2="$QCDIR/${s}_R2_fastqc.zip"

    if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
        echo "FASTQ files for sample $s do not exist. Skipping..."
        return
    fi

    if [[ -f "$output_r1" && -f "$output_r2" ]]; then
        echo "FastQC for sample $s already completed. Skipping..."
    else
        echo "Running FastQC for sample $s..."
        fastqc --quiet "$r1_file" --outdir="$QCDIR" -t 2 && echo "${s} read1 FastQC complete"
        fastqc --quiet "$r2_file" --outdir="$QCDIR" -t 2 && echo "${s} read2 FastQC complete"
    fi
}

export -f run_fastqc
export QCDIR READDIR

printf "%s\n" "${SAMPLES[@]}" | xargs -P 20 -I {} bash -c 'run_fastqc "$@"' _ {}

echo "Initial FastQC complete"
echo "Summarising FastQC results with MultiQC..."
multiqc "$QCDIR/" -o "$QCDIR" --force
echo "MultiQC report complete"


#######################################################################################
##                 2. Trimming and QC with Trim Galore (cutadapt + FastQC)           ##
#######################################################################################

TRIMDIR=$MAINDIR/trim_galore
mkdir -p "$TRIMDIR"

echo "Trimming reads..."

run_trim_galore() {
    local s="$1"
    local input_r1="$READDIR/${s}_R1.fastq.gz"
    local input_r2="$READDIR/${s}_R2.fastq.gz"
    local output_r1="$TRIMDIR/${s}_R1_val_1.fq.gz"
    local output_r2="$TRIMDIR/${s}_R2_val_2.fq.gz"

    if [[ ! -f "$input_r1" || ! -f "$input_r2" ]]; then
        echo "Input files for sample $s do not exist. Skipping..."
        return
    fi

    if [[ -f "$output_r1" && -f "$output_r2" ]]; then
        echo "Trim Galore for sample $s already completed. Skipping..."
    else
        echo "Trimming sample $s..."
        trim_galore --phred33 --paired --fastqc --cores 8 --output_dir "$TRIMDIR" \
            "$input_r1" "$input_r2"
    fi
}

export -f run_trim_galore
export TRIMDIR READDIR

printf "%s\n" "${SAMPLES[@]}" | xargs -P 5 -I {} bash -c 'run_trim_galore "$@"' _ {}
wait

echo "Trimming complete"
echo "Summarising post-trim FastQC results with MultiQC..."
multiqc "$TRIMDIR/" -o "$TRIMDIR" --force
echo "MultiQC report complete"


#######################################################################################
##                    3. Alignment to combined mm39 + dm6 genome                     ##
#######################################################################################

# Reads are aligned to a combined mm39/dm6 genome. dm6 chromosomes are prefixed with
# "dm6_" to allow clean separation after alignment. This spike-in strategy enables
# calibration of ChIP-seq signal across samples.

ALIGNDIR=$MAINDIR/bowtie2
mkdir -p "$ALIGNDIR"

echo "Aligning reads to combined mm39/dm6 genome (threads: $THREADS)..."

for s in "${SAMPLES[@]}"
do
    input_r1="$TRIMDIR/${s}_R1_val_1.fq.gz"
    input_r2="$TRIMDIR/${s}_R2_val_2.fq.gz"
    output_bam="$ALIGNDIR/${s}.bam"
    output_txt="$ALIGNDIR/${s}_bowtie2.txt"
    output_idxstats="$ALIGNDIR/${s}.idxstats"
    output_flagstat="$ALIGNDIR/${s}.flagstat"

    if [[ ! -f "$input_r1" || ! -f "$input_r2" ]]; then
        echo "Trimmed files for sample $s do not exist. Skipping..."
        continue
    fi

    if [[ -f "$output_bam" && -f "$output_txt" && -f "$output_idxstats" && -f "$output_flagstat" ]]; then
        echo "Alignment for sample $s already completed. Skipping..."
        continue
    fi

    echo "Aligning sample $s..."
    bowtie2 --very-sensitive-local --no-mixed --dovetail --phred33 \
        -I 10 -X 1000 \
        -p "$THREADS" \
        -x "$INDEX_PREFIX" \
        -1 "$input_r1" -2 "$input_r2" 2>"$output_txt" \
        | samtools sort -@ "$THREADS" -O bam -o "$output_bam"
    samtools index -@ "$THREADS" "$output_bam"
    samtools idxstats "$output_bam" > "$output_idxstats"
    samtools flagstat "$output_bam" > "$output_flagstat"
done
wait

echo "Alignment complete"
echo "Summarising alignment results with MultiQC..."
multiqc "$ALIGNDIR/" -o "$ALIGNDIR" --force
echo "MultiQC report complete"


#######################################################################################
##                    4. Separate mm39 and dm6 reads                                 ##
#######################################################################################

# dm6 reads are identified by the dm6_chr prefix added during index preparation.
# mm39 reads are extracted by negative selection (everything without the dm6_ prefix).

DM6DIR=$MAINDIR/bam_dm6
mkdir -p "$DM6DIR"

# --- 4a. Extract dm6 (spike-in) reads ---
echo "Extracting dm6 spike-in reads..."

for s in "${SAMPLES[@]}"
do
    input_bam="$ALIGNDIR/${s}.bam"
    output_bam="$DM6DIR/${s}_dm6.bam"

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: $input_bam not found. Skipping sample $s."
        continue
    fi

    if [[ -f "$output_bam" ]]; then
        echo "dm6 BAM for sample $s already exists. Skipping..."
        continue
    fi

    samtools idxstats "$input_bam" | cut -f 1 | grep dm6 \
        | xargs samtools view -b "$input_bam" > "$output_bam"
    samtools index "$output_bam"
done
wait

echo "dm6 reads extracted"

# --- 4b. Extract mm39 reads by negative selection ---
echo "Extracting mm39 reads..."

for s in "${SAMPLES[@]}"
do
    input_bam="$ALIGNDIR/${s}.bam"
    output_bam="$ALIGNDIR/${s}_mm39.bam"

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: $input_bam not found. Skipping sample $s."
        continue
    fi

    if [[ -f "$output_bam" ]]; then
        echo "mm39 BAM for sample $s already exists. Skipping..."
        continue
    fi

    samtools idxstats "$input_bam" | cut -f 1 | grep -v dm6 \
        | xargs samtools view -b "$input_bam" > "$output_bam"
    samtools index "$output_bam"
done
wait

echo "mm39 reads extracted"


#######################################################################################
##               5. Filtering alignment files with SAMtools and Picard               ##
#######################################################################################

FILTERDIR=$MAINDIR/cleanbams
mkdir -p "$FILTERDIR"

# --- 5a. Add read groups ---
# Alignment produces BAMs with no ReadGroup information; required by MarkDuplicates.

echo "Adding read groups with Picard..."

for s in "${SAMPLES[@]}"; do
    input_bam="$ALIGNDIR/${s}_mm39.bam"
    output_bam="$FILTERDIR/${s}_mm39_rg.bam"

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: $input_bam not found. Skipping sample $s."
        continue
    fi

    if [[ -f "$output_bam" ]]; then
        echo "Read groups for sample $s already added. Skipping..."
        continue
    fi

    picard AddOrReplaceReadGroups \
        -I "$input_bam" \
        -O "$output_bam" \
        -RGID "FLOWCELL1" \
        -RGLB "$s" \
        -RGPL "ILLUMINA" \
        -RGPU "run1" \
        -RGSM "$s"
done
wait

# --- 5b. Mark duplicates ---
# Note: for ChIP-seq, duplicates are marked but NOT removed (--REMOVE_DUPLICATES false)
# Duplicate rates are informative QC metrics for ChIP efficiency

echo "Marking duplicates with Picard..."

for s in "${SAMPLES[@]}"; do
    input_bam="$FILTERDIR/${s}_mm39_rg.bam"
    output_bam="$FILTERDIR/${s}_markdup.bam"

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: $input_bam not found. Skipping sample $s."
        continue
    fi

    if [[ -f "$output_bam" ]]; then
        echo "Duplicate marking for sample $s already completed. Skipping..."
        continue
    fi

    picard MarkDuplicates --QUIET true \
        --INPUT "$input_bam" \
        --OUTPUT "$output_bam" \
        --METRICS_FILE "$FILTERDIR/${s}.dup.metrics" \
        --REMOVE_DUPLICATES false \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT
done
wait

echo "Duplicate marking complete"
echo "Summarising deduplication results with MultiQC..."
multiqc "$FILTERDIR/" -o "$FILTERDIR" --force
echo "MultiQC report complete"


#######################################################################################
##                     6. Insert size metrics with Picard                            ##
#######################################################################################

INSERTDIR=$MAINDIR/insertmetrics
mkdir -p "$INSERTDIR"

echo "Collecting insert size metrics..."

for s in "${SAMPLES[@]}"
do
    input_bam="$FILTERDIR/${s}_markdup.bam"
    output_metrics="$INSERTDIR/${s}_insertMetrics.txt"
    output_hist="$INSERTDIR/${s}_insertSize.pdf"

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: $input_bam not found. Skipping sample $s."
        continue
    fi

    if [[ -f "$output_metrics" ]]; then
        echo "Insert size metrics for sample $s already exist. Skipping..."
        continue
    fi

    picard CollectInsertSizeMetrics --QUIET true \
        --INPUT "$input_bam" \
        --OUTPUT "$output_metrics" \
        --Histogram_FILE "$output_hist"
done
wait

echo "Insert size metrics complete"
echo "Pipeline complete. Clean BAMs are in: $FILTERDIR"
