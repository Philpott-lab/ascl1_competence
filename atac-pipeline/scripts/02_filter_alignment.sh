#!/bin/bash
set -euo pipefail

#######################################################################################
##                      ATAC-seq DATA PRE-PROCESSING & ALIGNMENT                     ##
#######################################################################################


#### OVERVIEW ####
# 1. FastQC checks the quality of the raw fastq sequencing reads
# 2. Trim Galore uses cutadapt and FastQC to trim adapters and check quality
# 3. Trimmed reads are aligned to the genome by bowtie2
# 4. Alignment files are filtered: mtDNA removal, read group addition, deduplication
# 5. Insert size metrics are collected using Picard


#### INPUT DATA ####
# Samples from mESCs with Ascl1-GR-HA
# 3 conditions: mESC, EpiLC, NE
# 3 timepoints per condition: 0h, 6h, 24h after dexamethasone induction of Ascl1
# Sample sheet: first column is the full sample name e.g. C11_EpiLC_HA_1
# FASTQ files expected as: <MAINDIR>/fastq/<sample>_r1.fq.gz and <sample>_r2.fq.gz
# Bowtie2 index expected at: <MAINDIR>/indexes/mm39 (run 01_prepare_index.sh first)


#### CONDA ENVIRONMENT ####
# conda create -n atacseq
# conda install -n atacseq \
#   cutadapt fastqc multiqc trim-galore bowtie2 samtools sambamba \
#   picard bedtools bamtools preseq pysam ucsc-bedgraphtobigwig \
#   deeptools homer subread
# See environment.yml for pinned versions.


#### USAGE ####
# bash 02_atac_preprocess_align.sh <main_directory> <sample_sheet> [threads]
#
# Arguments:
#   $1  MAINDIR     - path to the main project directory
#   $2  SAMPLESHEET - path to CSV sample sheet (first column = sample name, with header)
#   $3  THREADS     - number of threads for alignment and samtools (default: 8)
#
# Example:
#   bash 02_atac_preprocess_align.sh /path/to/project /path/to/samplesheet.csv 40


#### ARGUMENT HANDLING ####
if [[ $# -lt 2 ]]; then
    echo "Usage: bash 02_atac_preprocess_align.sh <main_directory> <sample_sheet> [threads]"
    exit 1
fi

MAINDIR=$1
SAMPLESHEET=$2
THREADS=${3:-8}

echo "Main directory : $MAINDIR"
echo "Sample sheet   : $SAMPLESHEET"
echo "Threads        : $THREADS"


#### CHECK BOWTIE2 INDEX EXISTS ####
INDEXDIR=$MAINDIR/indexes
if [[ ! -f "$INDEXDIR/mm39.1.bt2" ]]; then
    echo "[ERROR]: Bowtie2 index not found at $INDEXDIR/mm39."
    echo "         Please run 01_prepare_index.sh first."
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
    local r1_file="$READDIR/${s}_r1.fq.gz"
    local r2_file="$READDIR/${s}_r2.fq.gz"
    local output_r1="$QCDIR/${s}_r1_fastqc.zip"
    local output_r2="$QCDIR/${s}_r2_fastqc.zip"

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
    local input_r1="$READDIR/${s}_r1.fq.gz"
    local input_r2="$READDIR/${s}_r2.fq.gz"
    local output_r1="$TRIMDIR/${s}_r1_val_1.fq.gz"
    local output_r2="$TRIMDIR/${s}_r2_val_2.fq.gz"

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
##                           3. Alignment with Bowtie2                               ##
#######################################################################################

ALIGNDIR=$MAINDIR/bowtie2
mkdir -p "$ALIGNDIR"

echo "Aligning reads with Bowtie2 (threads: $THREADS)..."

for s in "${SAMPLES[@]}"
do
    input_r1="$TRIMDIR/${s}_r1_val_1.fq.gz"
    input_r2="$TRIMDIR/${s}_r2_val_2.fq.gz"
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

    echo "Aligning sample $s to mm39..."
    bowtie2 --very-sensitive-local --no-mixed --no-discordant -X 2000 \
        -p "$THREADS" -x "$INDEXDIR/mm39" \
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
##               4. Filtering alignment files with SAMtools and Picard               ##
#######################################################################################

# --- 4a. Calculate mtDNA content ---

output_file="$ALIGNDIR/mtDNA_content.txt"
echo "Sample mtDNA_Content(%)" > "$output_file"

echo "Calculating mtDNA content..."

for s in "${SAMPLES[@]}"
do
    bam_file="$ALIGNDIR/${s}.bam"

    if [[ ! -f "$bam_file" ]]; then
        echo "[ERROR]: BAM file $bam_file does not exist. Skipping."
        continue
    fi

    mtReads=$(samtools idxstats "$bam_file" | grep 'chrM' | cut -f 3)
    totalReads=$(samtools idxstats "$bam_file" | awk '{SUM += $3} END {print SUM}')

    if [[ -z "$mtReads" || -z "$totalReads" ]]; then
        echo "[ERROR]: Failed to calculate reads for $bam_file. Skipping."
        continue
    fi

    mtDNA_content=$(bc <<< "scale=2;100*$mtReads/$totalReads")
    echo "==> ${s} mtDNA content: ${mtDNA_content}%"
    echo "${s} ${mtDNA_content}%" >> "$output_file"
done

echo "mtDNA content written to $output_file"


# --- 4b. Remove mitochondrial reads ---

echo "Removing mitochondrial reads..."

remove_chrM() {
    local s="$1"
    local input_file="$ALIGNDIR/${s}.bam"
    local output_file="$ALIGNDIR/${s}.nochrM.bam"

    if [[ -f "$output_file" ]]; then
        echo "chrM-removed BAM for sample $s already exists. Skipping..."
        return 0
    fi

    if [[ -f "$input_file" ]]; then
        samtools view -h "$input_file" | grep -v chrM | samtools sort -O bam -o "$output_file"
        echo "chrM removal complete for sample $s"
    else
        echo "[ERROR]: $input_file not found. Skipping sample $s."
    fi
}

export -f remove_chrM
export ALIGNDIR

printf "%s\n" "${SAMPLES[@]}" | xargs -P 10 -I {} bash -c 'remove_chrM "$@"' _ {}
wait


# --- 4c. Add read groups with Picard ---
# Alignment produces BAMs with no ReadGroup information; required by MarkDuplicates.

FILTERDIR=$MAINDIR/cleanbams
mkdir -p "$FILTERDIR"

echo "Adding read groups with Picard..."

for s in "${SAMPLES[@]}"; do
    input_bam="$ALIGNDIR/${s}.nochrM.bam"
    output_bam="$FILTERDIR/${s}_rg.bam"

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
        -RGID "FLOWCELL1.LANE4" \
        -RGLB "$s" \
        -RGPL "ILLUMINA" \
        -RGPU "run1" \
        -RGSM "$s"
done
wait


# --- 4d. Mark and remove duplicates with Picard ---

echo "Marking and removing duplicates with Picard..."

for s in "${SAMPLES[@]}"; do
    input_bam="$FILTERDIR/${s}_rg.bam"
    output_bam="$FILTERDIR/${s}.nodup.bam"

    if [[ ! -f "$input_bam" ]]; then
        echo "[ERROR]: $input_bam not found. Skipping sample $s."
        continue
    fi

    if [[ -f "$output_bam" ]]; then
        echo "Deduplication for sample $s already completed. Skipping..."
        continue
    fi

    picard MarkDuplicates --QUIET true \
        --INPUT "$input_bam" \
        --OUTPUT "$output_bam" \
        --METRICS_FILE "$FILTERDIR/${s}.dup.metrics" \
        --REMOVE_DUPLICATES true \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT
done
wait

echo "Deduplication complete"
echo "Summarising deduplication results with MultiQC..."
multiqc "$FILTERDIR/" -o "$FILTERDIR" --force
echo "MultiQC report complete"


#######################################################################################
##                     5. Insert size metrics with Picard                            ##
#######################################################################################

INSERTDIR=$MAINDIR/insertmetrics
mkdir -p "$INSERTDIR"

echo "Collecting insert size metrics..."

for s in "${SAMPLES[@]}"
do
    input_bam="$FILTERDIR/${s}.nodup.bam"
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
