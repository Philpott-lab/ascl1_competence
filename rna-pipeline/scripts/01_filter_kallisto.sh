#!/bin/bash

MAINDIR=/data/20241016_LundieJ_IntegratedAnalysis/Analysis/data_rnaseq/preprocessing
INDEXDIR=$MAINDIR/index
READDIR=$MAINDIR/reads
OUTDIR=$MAINDIR/kallisto
QCDIR=$MAINDIR/initialfastqc
TRIMDIR=$MAINDIR/trimming

samples=(
    "exp1_gG5_mESC_t24_rep1" "exp1_gG5_mESC_t24_rep2" "exp1_gG5_mESC_t24_rep3"
    "exp1_gG5_EpiLC_t24_rep1" "exp1_gG5_EpiLC_t24_rep2" "exp1_gG5_EpiLC_t24_rep3"
    "exp1_gG5_NE_t24_rep1" "exp1_gG5_NE_t24_rep2" "exp1_gG5_NE_t24_rep3"
    "exp1_C11_mESC_t24_rep1" "exp1_C11_mESC_t24_rep2" "exp1_C11_mESC_t24_rep3"
    "exp1_C11_EpiLC_t24_rep1" "exp1_C11_EpiLC_t24_rep2" "exp1_C11_EpiLC_t24_rep3"
    "exp1_C11_NE_t24_rep1" "exp1_C11_NE_t24_rep2" "exp1_C11_NE_t24_rep3"
    "exp2_gG5_EpiLC_t0_rep1" "exp2_gG5_EpiLC_t0_rep2" "exp2_gG5_EpiLC_t0_rep3" "exp2_gG5_EpiLC_t0_rep4"
    "exp2_gG5_EpiLC_t24_rep1" "exp2_gG5_EpiLC_t24_rep2" "exp2_gG5_EpiLC_t24_rep3"
    "exp2_gG5_EpiLC_t6_rep1" "exp2_gG5_EpiLC_t6_rep2" "exp2_gG5_EpiLC_t6_rep3" "exp2_gG5_EpiLC_t6_rep4"
    "exp2_gG5_mESC_t0_rep1" "exp2_gG5_mESC_t0_rep2" "exp2_gG5_mESC_t0_rep3" "exp2_gG5_mESC_t0_rep4"
    "exp2_gG5_mESC_t24_rep1" "exp2_gG5_mESC_t24_rep2" "exp2_gG5_mESC_t24_rep3" "exp2_gG5_mESC_t24_rep4"
    "exp2_gG5_mESC_t6_rep1" "exp2_gG5_mESC_t6_rep2" "exp2_gG5_mESC_t6_rep3" "exp2_gG5_mESC_t6_rep4"
    "exp2_gG5_NE_t0_rep1" "exp2_gG5_NE_t0_rep2" "exp2_gG5_NE_t0_rep3" "exp2_gG5_NE_t0_rep4"
    "exp2_gG5_NE_t24_rep1" "exp2_gG5_NE_t24_rep2" "exp2_gG5_NE_t24_rep3" "exp2_gG5_NE_t24_rep4"
    "exp2_gG5_NE_t6_rep1" "exp2_gG5_NE_t6_rep2" "exp2_gG5_NE_t6_rep3" "exp2_gG5_NE_t6_rep4"
)

#######################################################################################
##                          1. Initial QC with FastQC                                ##
#######################################################################################

echo "Performing initial fastqc..."

# Define a function to run FastQC

run_fastqc() {
    local s="$1"
    local r1_file="$READDIR/${s}_r1.fq.gz"
    local r2_file="$READDIR/${s}_r2.fq.gz"
    local output_r1="$QCDIR/${s}_r1_fastqc.zip"
    local output_r2="$QCDIR/${s}_r2_fastqc.zip"

    # Check if the input files exist
    if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
        echo "Files for sample $s do not exist. Skipping..."
        return
    fi

    # Check if the output files already exist
    if [[ -f "$output_r1" && -f "$output_r2" ]]; then
        echo "FastQC for sample $s already completed. Skipping..."
    else
        echo "Running FastQC for sample $s..."
        fastqc --quiet "$r1_file" --outdir="$QCDIR" -t 2 && echo "${s} read1 FastQC Complete"
        fastqc --quiet "$r2_file" --outdir="$QCDIR" -t 2 && echo "${s} read2 FastQC Complete"
    fi
}

export -f run_fastqc
export QCDIR READDIR

# Run FastQC in parallel, with a maximum of 20 jobs at a time
printf "%s\n" "${samples[@]}" | xargs -P 20 -I {} bash -c 'run_fastqc "$@"' _ {}

echo "Initial FastQC Complete"
echo "Summarising FastQC results by MultiQC..."

multiqc "$QCDIR/" -o "$QCDIR" --force
echo "MultiQC report complete"




#######################################################################################
##                 2. Trimming and QC with trim_galore (cutadapt+FastQC)             ##
#######################################################################################

echo "Trimming reads..."

# Define the function for trim_galore
run_trim_galore() {
    local s="$1"
    local input_r1="$READDIR/${s}_r1.fq.gz"
    local input_r2="$READDIR/${s}_r2.fq.gz"
    local output_r1="$TRIMDIR/${s}_r1_val_1.fq.gz"
    local output_r2="$TRIMDIR/${s}_r2_val_2.fq.gz"
    
    # Check if the input files exist
    if [[ ! -f "$input_r1" || ! -f "$input_r2" ]]; then
        echo "Input files for sample $s do not exist. Skipping..."
        return
    fi

    # Check if the output files already exist
    if [[ -f "$output_r1" && -f "$output_r2" ]]; then
        echo "Files for sample $s already processed. Skipping..."
    else
        echo "Processing sample $s..."
        trim_galore --phred33 --paired --fastqc --cores 8 --output_dir "$TRIMDIR" \
        "$input_r1" \
        "$input_r2"
    fi
}

export -f run_trim_galore
export TRIMDIR READDIR

# Run Trim_Galore in parallel, with a maximum of 5 jobs at a time
printf "%s\n" "${samples[@]}" | xargs -P 8 -I {} bash -c 'run_trim_galore "$@"' _ {}
wait

echo "Trimming complete"

echo "Summarizing post-trim fastqc results by MultiQC..."

multiqc "$TRIMDIR/" -o "$TRIMDIR" --force

echo "multiqc report complete"



#######################################################################################
##                         3. Pseudoalignment with Kallisto                          ##
#######################################################################################

for samp in ${samples[@]}
do
    echo "Processing sample ${samp}"
    kallisto quant -i $INDEXDIR/index_mm39.idx -o $OUTDIR/${samp} -b 100 \
    --threads 16 \
    $TRIMDIR/${samp}_r1_val_1.fq.gz $TRIMDIR/${samp}_r2_val_2.fq.gz &> $OUTDIR/${samp}_kallisto.log
done
wait

multiqc "$OUTDIR/" -o "$OUTDIR" --force
