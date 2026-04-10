#!/bin/bash
# =============================================================================
# setup_ascl1_analysis.sh
# Philpott-lab/ascl1_competence — Reproducibility Setup Script
#
# Sets up the full pipeline environment on a new Linux machine.
# Assumes the following Data/ structure already exists (e.g. from GEO download):
#
#   <data-dir>/
#     atac-seq/    ← ATAC-seq FASTQ files
#     chip-seq/    ← ChIP-seq FASTQ files
#     rna-seq/     ← RNA-seq FASTQ files
#
# What this script does:
#   1. Clones the GitHub repo alongside Data/
#   2. Creates pipeline output subdirectories inside each data folder
#   3. Creates fastq/ symlinks so pipeline scripts find the FASTQs
#   4. Generates config.yaml — single source of truth for all parameters
#   5. Installs conda environments (atacseq, chipseq) and yq (YAML parser)
#   6. Installs required R/Bioconductor packages (including yaml)
#   7. Builds the Kallisto mm39 transcriptome index
#   8. Creates a patched RNA-seq script with correct paths
#   9. Generates run_pipeline.sh — reads all parameters from config.yaml
#
# Usage:
#   bash setup_ascl1_analysis.sh \
#       --install-dir   /path/to/install \
#       --data-dir      /path/to/Data \
#       [--threads      16] \
#       [--skip-conda]
#
#   --install-dir   Where to clone the repo. The repo lands at
#                   <install-dir>/ascl1_competence/ and run_pipeline.sh
#                   is written to <install-dir>/.
#   --data-dir      Path to your Data/ folder containing atac-seq/,
#                   chip-seq/, and rna-seq/ subdirectories.
#   --threads       CPU threads (default: 8)
#   --skip-conda    Skip conda env creation (use if envs already exist)
# =============================================================================

set -euo pipefail

# ─── Defaults ─────────────────────────────────────────────────────────────────
INSTALL_DIR=""
DATA_DIR=""
THREADS=8
SKIP_CONDA=false
REPO_URL="https://github.com/Philpott-lab/ascl1_competence.git"

# ─── Parse arguments ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-dir)  INSTALL_DIR="$2"; shift 2 ;;
        --data-dir)     DATA_DIR="$2";    shift 2 ;;
        --threads)      THREADS="$2";     shift 2 ;;
        --skip-conda)   SKIP_CONDA=true;  shift ;;
        *) echo "[ERROR]: Unknown argument: $1"; exit 1 ;;
    esac
done

# ─── Validate ─────────────────────────────────────────────────────────────────
missing=false
if [[ -z "$INSTALL_DIR" ]]; then echo "[ERROR]: --install-dir is required"; missing=true; fi
if [[ -z "$DATA_DIR" ]];    then echo "[ERROR]: --data-dir is required";    missing=true; fi
if [[ "$missing" == true ]]; then
    echo ""
    echo "Usage: bash setup_ascl1_analysis.sh \\"
    echo "    --install-dir  /path/to/install \\"
    echo "    --data-dir     /path/to/Data \\"
    echo "    [--threads 8] [--skip-conda]"
    exit 1
fi

INSTALL_DIR=$(realpath "$INSTALL_DIR")
DATA_DIR=$(realpath "$DATA_DIR")

# Check expected subdirectories exist in Data/
for subdir in atac-seq chip-seq rna-seq; do
    if [[ ! -d "$DATA_DIR/$subdir" ]]; then
        echo "[ERROR]: Expected directory not found: $DATA_DIR/$subdir"
        echo "         Data/ should contain: atac-seq/, chip-seq/, rna-seq/"
        exit 1
    fi
done

REPO_DIR="$INSTALL_DIR/ascl1_competence"
FIGURES_DIR="$INSTALL_DIR/Figures"
CONFIG_FILE="$REPO_DIR/config.yaml"

ATAC_DIR="$DATA_DIR/atac-seq"
CHIP_DIR="$DATA_DIR/chip-seq"
RNA_DIR="$DATA_DIR/rna-seq"

echo "============================================================"
echo " ascl1_competence analysis — environment setup"
echo "============================================================"
echo "  Install directory : $INSTALL_DIR"
echo "  Repo directory    : $REPO_DIR"
echo "  Data directory    : $DATA_DIR"
echo "  Config file       : $CONFIG_FILE"
echo "  ATAC data         : $ATAC_DIR"
echo "  ChIP data         : $CHIP_DIR"
echo "  RNA data          : $RNA_DIR"
echo "  Threads           : $THREADS"
echo "  Skip conda        : $SKIP_CONDA"
echo "============================================================"
echo ""

mkdir -p "$INSTALL_DIR" "$FIGURES_DIR"


# =============================================================================
# STEP 1: Clone the repository
# =============================================================================

echo "──────────────────────────────────────────────────────"
echo "STEP 1: Cloning repository"
echo "──────────────────────────────────────────────────────"

if [[ -d "$REPO_DIR/.git" ]]; then
    echo "Repository already exists. Pulling latest..."
    git -C "$REPO_DIR" pull
else
    git clone "$REPO_URL" "$REPO_DIR"
fi

echo "Repository ready at: $REPO_DIR"


# =============================================================================
# STEP 2: Create pipeline output directories and fastq/ symlinks
# =============================================================================
# The pipeline scripts expect FASTQs at $MAINDIR/fastq/<sample>.fq.gz.
# Since FASTQs live directly in atac-seq/, chip-seq/, and rna-seq/, we create
# a fastq/ symlink inside each pointing to the parent directory, so the
# pipeline scripts can find files without modification.
#
# Similarly the R scripts use here("../Data/...") relative to the repo root,
# so Data/ must sit alongside the repo (which it does if --install-dir and
# --data-dir share the same parent, or Data/ is symlinked there).

echo ""
echo "──────────────────────────────────────────────────────"
echo "STEP 2: Creating pipeline directories and symlinks"
echo "──────────────────────────────────────────────────────"

# Create output subdirectories
mkdir -p "$ATAC_DIR/indexes"
mkdir -p "$CHIP_DIR/indexes"
mkdir -p "$RNA_DIR/index"
for subexp in ascl1 hdaci phox2b; do
    mkdir -p "$RNA_DIR/$subexp/kallisto"
    mkdir -p "$RNA_DIR/$subexp/initialfastqc"
    mkdir -p "$RNA_DIR/$subexp/trimming"
done

# fastq/ symlinks — point back to the data directory itself so that
# $MAINDIR/fastq/<sample>.fq.gz resolves to $MAINDIR/<sample>.fq.gz
for dir in "$ATAC_DIR" "$CHIP_DIR"; do
    if [[ ! -L "$dir/fastq" && ! -d "$dir/fastq" ]]; then
        ln -s "$dir" "$dir/fastq"
        echo "Created fastq symlink: $dir/fastq → $dir"
    else
        echo "fastq symlink already exists in $dir. Skipping."
    fi
done

# If Data/ is not already a sibling of the repo, symlink it there
# so that here("../Data/...") in R resolves correctly.
EXPECTED_DATA="$INSTALL_DIR/Data"
if [[ "$DATA_DIR" != "$EXPECTED_DATA" ]]; then
    if [[ ! -L "$EXPECTED_DATA" && ! -d "$EXPECTED_DATA" ]]; then
        ln -s "$DATA_DIR" "$EXPECTED_DATA"
        echo "Symlinked Data/ into install dir: $EXPECTED_DATA → $DATA_DIR"
    else
        echo "Data/ symlink already exists at $EXPECTED_DATA. Skipping."
    fi
fi

# Copy DiffBind sample sheets into their expected Data/ locations
# R scripts read: here("../Data/atac-seq/atacseq_diffbindSamples.csv") etc.
cp "$REPO_DIR/atac-pipeline/atacseq_diffbindSamples.csv" \
   "$ATAC_DIR/atacseq_diffbindSamples.csv"
cp "$REPO_DIR/chip-pipeline/ChIPseq_ASCL1_samples.csv" \
   "$CHIP_DIR/ChIPseq_ASCL1_samples.csv"

echo "Pipeline directories and symlinks ready."


# =============================================================================
# STEP 3: Generate config.yaml
# =============================================================================

echo ""
echo "──────────────────────────────────────────────────────"
echo "STEP 3: Generating config.yaml"
echo "──────────────────────────────────────────────────────"

if [[ -f "$CONFIG_FILE" ]]; then
    echo "config.yaml already exists — preserving existing settings."
    echo "  Delete $CONFIG_FILE and re-run setup to regenerate."
else
    cat > "$CONFIG_FILE" << YAML
# =============================================================================
# config.yaml — ascl1_competence analysis global parameters
#
# Single source of truth for all pipeline parameters.
# Read by run_pipeline.sh (via yq) and R scripts (via the yaml package).
# Commit this file to git to document the parameters used for each run.
# =============================================================================

# ── Directory paths ──────────────────────────────────────────────────────────
paths:
  install_dir:  "$INSTALL_DIR"      # project root
  repo_dir:     "$REPO_DIR"         # git repository
  data_dir:     "$DATA_DIR"         # Data/ folder (atac-seq/, chip-seq/, rna-seq/)
  figures_dir:  "$FIGURES_DIR"      # output figures

# ── Compute resources ────────────────────────────────────────────────────────
compute:
  threads:       $THREADS   # CPUs for alignment, bamCoverage, etc.
  parallel_jobs: 4          # samples processed in parallel (peak calling)

# ── Genome parameters ────────────────────────────────────────────────────────
genome:
  mm39_size:       2468088461   # effective genome size for mm39
  ensembl_release: 113          # Ensembl release for Kallisto index

# ── ATAC-seq parameters ──────────────────────────────────────────────────────
atac:
  diffbind_min_overlap: 3    # min replicates per condition for consensus peaks
  diffbind_summits:     75   # ±bp extension around peak summit for counting
  bigwig_bin_size:      10   # bigwig bin size in bp

# ── ChIP-seq parameters ──────────────────────────────────────────────────────
chip:
  bigwig_bin_size: 10

# ── RNA-seq parameters ───────────────────────────────────────────────────────
rna:
  kallisto_bootstraps: 100
  # phox2b is single-end; estimate fragment length from FastQC/Picard before running
  phox2b_fragment_length: 200
  phox2b_fragment_sd:     20
YAML

    echo "config.yaml written to: $CONFIG_FILE"
fi


# =============================================================================
# STEP 4: Install conda environments and yq
# =============================================================================

echo ""
echo "──────────────────────────────────────────────────────"
echo "STEP 4: Installing conda environments and yq"
echo "──────────────────────────────────────────────────────"

if [[ "$SKIP_CONDA" == true ]]; then
    echo "Skipping (--skip-conda). Ensure atacseq/chipseq envs and yq are available."
else
    if ! command -v conda &> /dev/null; then
        echo "[ERROR]: conda not found. Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi

    if command -v yq &> /dev/null; then
        echo "yq already installed: $(yq --version)"
    else
        echo "Installing yq..."
        conda install -y -c conda-forge yq
    fi

    for env_name in atacseq chipseq; do
        [[ "$env_name" == "atacseq" ]] && YML="$REPO_DIR/atac-pipeline/env_atacseq.yml" \
                                       || YML="$REPO_DIR/chip-pipeline/env_chipseq.yml"
        YML_CLEAN="/tmp/${env_name}_clean.yml"
        grep -v '^prefix:' "$YML" > "$YML_CLEAN"   # strip machine-specific prefix line

        if conda env list | grep -q "^${env_name} "; then
            echo "Conda env '${env_name}' already exists. Skipping."
        else
            echo "Creating conda env: ${env_name} (may take several minutes)..."
            conda env create --name "$env_name" --file "$YML_CLEAN"
        fi
    done
fi


# =============================================================================
# STEP 5: Install R packages
# =============================================================================

echo ""
echo "──────────────────────────────────────────────────────"
echo "STEP 5: Installing R packages"
echo "──────────────────────────────────────────────────────"

if ! command -v Rscript &> /dev/null; then
    echo "[WARNING]: Rscript not found. Install R >=4.4, then install:"
    echo "  CRAN:        here, tidyverse, ggplot2, cowplot, yaml"
    echo "  Bioconductor: DiffBind, AnnotationHub, GenomicRanges, DESeq2,"
    echo "                ChIPseeker, TxDb.Mmusculus.UCSC.mm39.knownGene,"
    echo "                org.Mm.eg.db, clusterProfiler, PCAtools, edgeR, tximport"
else
    Rscript - <<'REOF'
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

cran_pkgs <- c("here", "tidyverse", "ggplot2", "cowplot", "yaml", "eulerr")
bioc_pkgs <- c("DiffBind", "AnnotationHub", "GenomicRanges", "DESeq2",
               "ChIPseeker", "TxDb.Mmusculus.UCSC.mm39.knownGene",
               "org.Mm.eg.db", "clusterProfiler", "PCAtools", "edgeR", "tximport",
               "BSgenome.Mmusculus.UCSC.mm39", "rtracklayer", "GenomicInteractions")

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing: ", pkg); install.packages(pkg, repos = "https://cloud.r-project.org")
    } else message("Already installed: ", pkg)
}
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing: ", pkg); BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else message("Already installed: ", pkg)
}
message("R packages ready.")
REOF
fi


# =============================================================================
# STEP 6: Build Kallisto mm39 index
# =============================================================================

echo ""
echo "──────────────────────────────────────────────────────"
echo "STEP 6: Building Kallisto mm39 transcriptome index"
echo "──────────────────────────────────────────────────────"

KALLISTO_INDEX="$RNA_DIR/index/index_mm39.idx"
MM39_CDNA="$RNA_DIR/index/Mus_musculus.GRCm39.cdna.all.fa.gz"

if [[ -f "$KALLISTO_INDEX" ]]; then
    echo "Kallisto index already exists. Skipping."
elif ! command -v kallisto &> /dev/null; then
    echo "[WARNING]: kallisto not found. Build manually after installing:"
    echo "  wget -P $RNA_DIR/index/ https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
    echo "  kallisto index -i $KALLISTO_INDEX $MM39_CDNA"
else
    [[ ! -f "$MM39_CDNA" ]] && wget -q -P "$RNA_DIR/index/" \
        "https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
    kallisto index -i "$KALLISTO_INDEX" "$MM39_CDNA"
    echo "Index built: $KALLISTO_INDEX"
fi


# =============================================================================
# STEP 7: Create patched RNA-seq script (reads paths from config.yaml)
# =============================================================================

echo ""
echo "──────────────────────────────────────────────────────"
echo "STEP 7: Creating RNA-seq scripts (ascl1, hdaci, phox2b)"
echo "──────────────────────────────────────────────────────"

# ── 7a: ascl1 experiment (exp1 + exp2, paired-end) ───────────────────────────
# Patched from 01_filter_kallisto.sh with paths redirected to rna-seq/ascl1/.

ASCL1_RNA="$INSTALL_DIR/run_rnaseq_ascl1.sh"

cat > "$ASCL1_RNA" << 'ASCL1HEADER'
#!/bin/bash
set -euo pipefail
# ascl1 experiment RNA-seq Kallisto script (exp1 + exp2, paired-end).
# Paths from config.yaml. Generated by setup_ascl1_analysis.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/ascl1_competence/config.yaml"

MAINDIR="$(yq '.paths.data_dir' "$CONFIG")/rna-seq/ascl1"
INDEXDIR="$(yq '.paths.data_dir' "$CONFIG")/rna-seq/index"
READDIR="$MAINDIR"
OUTDIR="$MAINDIR/kallisto"
QCDIR="$MAINDIR/initialfastqc"
TRIMDIR="$MAINDIR/trimming"
KALLISTO_BOOTSTRAPS="$(yq '.rna.kallisto_bootstraps' "$CONFIG")"
THREADS="$(yq '.compute.threads' "$CONFIG")"

mkdir -p "$OUTDIR" "$QCDIR" "$TRIMDIR"

ASCL1HEADER

tail -n +10 "$REPO_DIR/rna-pipeline/scripts/01_filter_kallisto.sh" \
    | grep -v '^MAINDIR=' | grep -v '^INDEXDIR=' | grep -v '^READDIR=' \
    | grep -v '^OUTDIR='  | grep -v '^QCDIR='    | grep -v '^TRIMDIR=' \
    | sed 's/-b 100/-b $KALLISTO_BOOTSTRAPS/g' \
    | sed 's/--threads 16/--threads $THREADS/g' \
    >> "$ASCL1_RNA"

chmod +x "$ASCL1_RNA"
echo "ascl1 RNA-seq script: $ASCL1_RNA"


# ── 7b: hdaci experiment (exp3, paired-end, dynamic sample discovery) ────────
# Sample list is discovered at script-run time so that the hdaci directory can
# be populated after setup runs (file copying may still be in progress).

HDACI_RNA="$INSTALL_DIR/run_rnaseq_hdaci.sh"

cat > "$HDACI_RNA" << 'HDACIEOF'
#!/bin/bash
set -euo pipefail
# hdaci experiment RNA-seq Kallisto script (exp3, paired-end).
# Samples are discovered dynamically from files present in rna-seq/hdaci/.
# Paths from config.yaml. Generated by setup_ascl1_analysis.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/ascl1_competence/config.yaml"

MAINDIR="$(yq '.paths.data_dir' "$CONFIG")/rna-seq/hdaci"
INDEXDIR="$(yq '.paths.data_dir' "$CONFIG")/rna-seq/index"
READDIR="$MAINDIR"
OUTDIR="$MAINDIR/kallisto"
QCDIR="$MAINDIR/initialfastqc"
TRIMDIR="$MAINDIR/trimming"
KALLISTO_BOOTSTRAPS="$(yq '.rna.kallisto_bootstraps' "$CONFIG")"
THREADS="$(yq '.compute.threads' "$CONFIG")"

mkdir -p "$OUTDIR" "$QCDIR" "$TRIMDIR"

# Discover samples from *_r1.fq.gz files present in the hdaci directory
mapfile -t samples < <(
    find "$READDIR" -maxdepth 1 -name '*_r1.fq.gz' \
    | sed 's|.*/||; s/_r1\.fq\.gz$//' | sort
)

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "[ERROR]: No *_r1.fq.gz files found in $READDIR"
    exit 1
fi
echo "Found ${#samples[@]} hdaci samples."

#######################################################################################
##                          1. Initial QC with FastQC                                ##
#######################################################################################

echo "Performing initial fastqc..."

run_fastqc() {
    local s="$1"
    local r1_file="$READDIR/${s}_r1.fq.gz"
    local r2_file="$READDIR/${s}_r2.fq.gz"
    local output_r1="$QCDIR/${s}_r1_fastqc.zip"
    local output_r2="$QCDIR/${s}_r2_fastqc.zip"
    if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
        echo "Files for sample $s do not exist. Skipping..."
        return
    fi
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
printf "%s\n" "${samples[@]}" | xargs -P 20 -I {} bash -c 'run_fastqc "$@"' _ {}
echo "Initial FastQC Complete"
multiqc "$QCDIR/" -o "$QCDIR" --force
echo "MultiQC report complete"

#######################################################################################
##                 2. Trimming and QC with trim_galore (cutadapt+FastQC)             ##
#######################################################################################

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
        echo "Files for sample $s already processed. Skipping..."
    else
        echo "Processing sample $s..."
        trim_galore --phred33 --paired --fastqc --cores 8 --output_dir "$TRIMDIR" \
            "$input_r1" "$input_r2"
    fi
}
export -f run_trim_galore
export TRIMDIR READDIR
printf "%s\n" "${samples[@]}" | xargs -P 8 -I {} bash -c 'run_trim_galore "$@"' _ {}
wait
echo "Trimming complete"
multiqc "$TRIMDIR/" -o "$TRIMDIR" --force
echo "multiqc report complete"

#######################################################################################
##                         3. Pseudoalignment with Kallisto                          ##
#######################################################################################

for samp in "${samples[@]}"; do
    echo "Processing sample ${samp}"
    kallisto quant -i "$INDEXDIR/index_mm39.idx" -o "$OUTDIR/${samp}" \
        -b "$KALLISTO_BOOTSTRAPS" --threads "$THREADS" \
        "$TRIMDIR/${samp}_r1_val_1.fq.gz" "$TRIMDIR/${samp}_r2_val_2.fq.gz" \
        &> "$OUTDIR/${samp}_kallisto.log"
done
wait
multiqc "$OUTDIR/" -o "$OUTDIR" --force
HDACIEOF

chmod +x "$HDACI_RNA"
echo "hdaci RNA-seq script: $HDACI_RNA"


# ── 7c: phox2b experiment (single-end, .fastq.gz) ────────────────────────────
# These are single-end reads; kallisto requires --single -l <len> -s <sd>.
# Fragment length/SD are set in config.yaml under rna.phox2b_fragment_length/sd.
# Estimate these values from FastQC before running (see config.yaml comments).

PHOX2B_RNA="$INSTALL_DIR/run_rnaseq_phox2b.sh"

cat > "$PHOX2B_RNA" << 'PHOX2BEOF'
#!/bin/bash
set -euo pipefail
# phox2b experiment RNA-seq Kallisto script (single-end, .fastq.gz).
# Paths and fragment parameters from config.yaml.
# Generated by setup_ascl1_analysis.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/ascl1_competence/config.yaml"

MAINDIR="$(yq '.paths.data_dir' "$CONFIG")/rna-seq/phox2b"
INDEXDIR="$(yq '.paths.data_dir' "$CONFIG")/rna-seq/index"
READDIR="$MAINDIR"
OUTDIR="$MAINDIR/kallisto"
QCDIR="$MAINDIR/initialfastqc"
TRIMDIR="$MAINDIR/trimming"
KALLISTO_BOOTSTRAPS="$(yq '.rna.kallisto_bootstraps' "$CONFIG")"
FRAG_LEN="$(yq '.rna.phox2b_fragment_length' "$CONFIG")"
FRAG_SD="$(yq '.rna.phox2b_fragment_sd' "$CONFIG")"
THREADS="$(yq '.compute.threads' "$CONFIG")"

mkdir -p "$OUTDIR" "$QCDIR" "$TRIMDIR"

samples=(
    "Z36K4W_1_Phox2b_mESC_NoAscl1_NoDox_Rep_1"
    "Z36K4W_2_Phox2b_mESC_NoAscl1_NoDox_Rep_2"
    "Z36K4W_3_Phox2b_mESC_NoAscl1_NoDox_Rep_3"
    "Z36K4W_4_Phox2b_mESC_Ascl1_NoDox_Rep_1"
    "Z36K4W_5_Phox2b_mESC_Ascl1_NoDox_Rep_2"
    "Z36K4W_6_Phox2b_mESC_Ascl1_NoDox_Rep_3"
    "Z36K4W_7_Phox2b_mESC_NoAscl1_Dox_Rep_1"
    "Z36K4W_8_Phox2b_mESC_NoAscl1_Dox_Rep_2"
    "Z36K4W_9_Phox2b_mESC_NoAscl1_Dox_Rep_3"
    "Z36K4W_10_Phox2b_mESC_Ascl1_Dox_Rep_1"
    "Z36K4W_11_Phox2b_mESC_Ascl1_Dox_Rep_2"
    "Z36K4W_12_Phox2b_mESC_Ascl1_Dox_Rep_3"
    "Z36K4W_13_Phox2b_EpiLC_NoAscl1_NoDox_Rep_1"
    "Z36K4W_14_Phox2b_EpiLC_NoAscl1_NoDox_Rep_2"
    "Z36K4W_15_Phox2b_EpiLC_NoAscl1_NoDox_Rep_3"
    "Z36K4W_16_Phox2b_EpiLC_Ascl1_NoDox_Rep_1"
    "Z36K4W_17_Phox2b_EpiLC_Ascl1_NoDox_Rep_2"
    "Z36K4W_18_Phox2b_EpiLC_Ascl1_NoDox_Rep_3"
    "Z36K4W_19_Phox2b_EpiLC_NoAscl1_Dox_Rep_1"
    "Z36K4W_20_Phox2b_EpiLC_NoAscl1_Dox_Rep_2"
    "Z36K4W_21_Phox2b_EpiLC_NoAscl1_Dox_Rep_3"
    "Z36K4W_22_Phox2b_EpiLC_Ascl1_Dox_Rep_1"
    "Z36K4W_23_Phox2b_EpiLC_Ascl1_Dox_Rep_2"
    "Z36K4W_24_Phox2b_EpiLC_Ascl1_Dox_Rep_3"
    "Z36K4W_25_Phox2b_NE_NoAscl1_NoDox_Rep_1"
    "Z36K4W_26_Phox2b_NE_NoAscl1_NoDox_Rep_2"
    "Z36K4W_27_Phox2b_NE_NoAscl1_NoDox_Rep_3"
    "Z36K4W_28_Phox2b_NE_Ascl1_NoDox_Rep_1"
    "Z36K4W_29_Phox2b_NE_Ascl1_NoDox_Rep_2"
    "Z36K4W_30_Phox2b_NE_Ascl1_NoDox_Rep_3"
    "Z36K4W_31_Phox2b_NE_NoAscl1_Dox_Rep_1"
    "Z36K4W_32_Phox2b_NE_NoAscl1_Dox_Rep_2"
    "Z36K4W_33_Phox2b_NE_NoAscl1_Dox_Rep_3"
    "Z36K4W_34_Phox2b_NE_Ascl1_Dox_Rep_1"
    "Z36K4W_35_Phox2b_NE_Ascl1_Dox_Rep_2"
    "Z36K4W_36_Phox2b_NE_Ascl1_Dox_Rep_3"
)

#######################################################################################
##                          1. Initial QC with FastQC                                ##
#######################################################################################

echo "Performing initial fastqc..."
for s in "${samples[@]}"; do
    input="$READDIR/${s}.fastq.gz"
    output="$QCDIR/${s}_fastqc.zip"
    if [[ ! -f "$input" ]]; then
        echo "File for sample $s does not exist. Skipping..."
    elif [[ -f "$output" ]]; then
        echo "FastQC for sample $s already completed. Skipping..."
    else
        echo "Running FastQC for sample $s..."
        fastqc --quiet "$input" --outdir="$QCDIR" -t 2
    fi
done
multiqc "$QCDIR/" -o "$QCDIR" --force
echo "Initial FastQC Complete"

#######################################################################################
##                 2. Trimming and QC with trim_galore (single-end)                  ##
#######################################################################################

echo "Trimming reads..."
for s in "${samples[@]}"; do
    input="$READDIR/${s}.fastq.gz"
    output="$TRIMDIR/${s}_trimmed.fastq.gz"
    if [[ ! -f "$input" ]]; then
        echo "Input file for sample $s does not exist. Skipping..."
    elif [[ -f "$output" ]]; then
        echo "Files for sample $s already processed. Skipping..."
    else
        echo "Processing sample $s..."
        trim_galore --phred33 --fastqc --cores 8 --output_dir "$TRIMDIR" "$input"
    fi
done
multiqc "$TRIMDIR/" -o "$TRIMDIR" --force
echo "Trimming complete"

#######################################################################################
##               3. Pseudoalignment with Kallisto (single-end)                       ##
##  Fragment length/SD: set rna.phox2b_fragment_length/sd in config.yaml            ##
#######################################################################################

for samp in "${samples[@]}"; do
    echo "Processing sample ${samp}"
    kallisto quant -i "$INDEXDIR/index_mm39.idx" -o "$OUTDIR/${samp}" \
        -b "$KALLISTO_BOOTSTRAPS" --single -l "$FRAG_LEN" -s "$FRAG_SD" \
        --threads "$THREADS" \
        "$TRIMDIR/${samp}_trimmed.fastq.gz" \
        &> "$OUTDIR/${samp}_kallisto.log"
done
wait
multiqc "$OUTDIR/" -o "$OUTDIR" --force
PHOX2BEOF

chmod +x "$PHOX2B_RNA"
echo "phox2b RNA-seq script: $PHOX2B_RNA"


# =============================================================================
# STEP 8: Generate run_pipeline.sh
# =============================================================================

echo ""
echo "──────────────────────────────────────────────────────"
echo "STEP 8: Generating run_pipeline.sh"
echo "──────────────────────────────────────────────────────"

RUN_SCRIPT="$INSTALL_DIR/run_pipeline.sh"

cat > "$RUN_SCRIPT" << 'RUNSCRIPT'
#!/bin/bash
# =============================================================================
# run_pipeline.sh — ascl1_competence master pipeline runner
#
# All parameters are read from ascl1_competence/config.yaml.
# Edit config.yaml to change any setting — do not edit this file directly.
#
# Usage:
#   bash run_pipeline.sh              # run all three pipelines
#   bash run_pipeline.sh --atac-only
#   bash run_pipeline.sh --chip-only
#   bash run_pipeline.sh --rna-only
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/ascl1_competence/config.yaml"

if [[ ! -f "$CONFIG" ]]; then
    echo "[ERROR]: config.yaml not found at $CONFIG"; exit 1
fi
if ! command -v yq &> /dev/null; then
    echo "[ERROR]: yq not found. Install: conda install -c conda-forge yq"; exit 1
fi

# ── Read parameters from config.yaml ─────────────────────────────────────────
INSTALL_DIR="$(yq  '.paths.install_dir'         "$CONFIG")"
REPO_DIR="$(yq     '.paths.repo_dir'            "$CONFIG")"
DATA_DIR="$(yq     '.paths.data_dir'            "$CONFIG")"

THREADS="$(yq      '.compute.threads'           "$CONFIG")"
PARALLEL_JOBS="$(yq '.compute.parallel_jobs'    "$CONFIG")"
MM39_SIZE="$(yq    '.genome.mm39_size'          "$CONFIG")"
ATAC_BIGWIG_BIN="$(yq '.atac.bigwig_bin_size'   "$CONFIG")"
CHIP_BIGWIG_BIN="$(yq '.chip.bigwig_bin_size'   "$CONFIG")"

ATAC_DIR="$DATA_DIR/atac-seq"
CHIP_DIR="$DATA_DIR/chip-seq"
ATAC_SCRIPTS="$REPO_DIR/atac-pipeline/scripts"
CHIP_SCRIPTS="$REPO_DIR/chip-pipeline/scripts"

# ── Flags ─────────────────────────────────────────────────────────────────────
RUN_ATAC=true; RUN_CHIP=true; RUN_RNA=true
for arg in "$@"; do
    case "$arg" in
        --atac-only) RUN_CHIP=false; RUN_RNA=false ;;
        --chip-only) RUN_ATAC=false; RUN_RNA=false ;;
        --rna-only)  RUN_ATAC=false; RUN_CHIP=false ;;
    esac
done

log() { echo "$(date '+%Y-%m-%d %H:%M:%S')  $1"; }

log "Config: $CONFIG"
log "  Threads: $THREADS  |  Parallel jobs: $PARALLEL_JOBS  |  mm39: $MM39_SIZE"


# =============================================================================
# ATAC-SEQ PIPELINE
# =============================================================================
if [[ "$RUN_ATAC" == true ]]; then

log "══════════════════════════ ATAC-SEQ ══════════════════════════"
eval "$(conda shell.bash hook)" && conda activate atacseq

log "ATAC 1/4: Building Bowtie2 mm39 index..."
bash "$ATAC_SCRIPTS/01_prepare_index.sh" "$ATAC_DIR" "$THREADS"

log "ATAC 2/4: Preprocessing and aligning reads..."
bash "$ATAC_SCRIPTS/02_filter_alignment.sh" \
    "$ATAC_DIR" \
    "$REPO_DIR/atac-pipeline/atacseq_bowtieSamples.csv" \
    "$THREADS"

# bam/ symlink: pipeline outputs to cleanbams/; DiffBind CSV references bam/
[[ ! -L "$ATAC_DIR/bam" ]] \
    && ln -s "$ATAC_DIR/cleanbams" "$ATAC_DIR/bam" \
    && log "Created symlink: atac-seq/bam → atac-seq/cleanbams"

log "ATAC 3/4: Calling peaks (MACS3, genome: $MM39_SIZE)..."
bash "$ATAC_SCRIPTS/03_callpeaks.sh" \
    "$ATAC_DIR" \
    "$REPO_DIR/atac-pipeline/atacseq_bowtieSamples.csv" \
    "$PARALLEL_JOBS"

# peaks/ symlink: pipeline outputs to macs3/; DiffBind CSV references peaks/
[[ ! -L "$ATAC_DIR/peaks" ]] \
    && ln -s "$ATAC_DIR/macs3" "$ATAC_DIR/peaks" \
    && log "Created symlink: atac-seq/peaks → atac-seq/macs3"

log "ATAC 4a/4: DiffBind TMM normalisation (R)..."
cd "$REPO_DIR" && Rscript "$ATAC_SCRIPTS/04a_diffbind_normalise.R" && cd "$INSTALL_DIR"

log "ATAC 4b/4: Generating TMM bigwigs (bin: ${ATAC_BIGWIG_BIN}bp)..."
bash "$ATAC_SCRIPTS/04b_bam_coverage.sh" \
    "$ATAC_DIR/cleanbams" \
    "$DATA_DIR/ATAC_TMM_scalefactors.csv" \
    "$THREADS" "$ATAC_BIGWIG_BIN"

log "ATAC pipeline complete."
fi


# =============================================================================
# CHIP-SEQ PIPELINE
# =============================================================================
if [[ "$RUN_CHIP" == true ]]; then

log "══════════════════════════ CHIP-SEQ ══════════════════════════"
eval "$(conda shell.bash hook)" && conda activate chipseq

log "ChIP 1/5: Building mm39+dm6 Bowtie2 index..."
bash "$CHIP_SCRIPTS/01_prepare_index.sh" "$CHIP_DIR" "$THREADS"

log "ChIP 2/5: Preprocessing and aligning reads..."
bash "$CHIP_SCRIPTS/02_preprocess_align.sh" \
    "$CHIP_DIR" \
    "$REPO_DIR/chip-pipeline/ChIPseq_align_samples.csv" \
    "$THREADS"

[[ ! -L "$CHIP_DIR/bam" ]] \
    && ln -s "$CHIP_DIR/cleanbams" "$CHIP_DIR/bam" \
    && log "Created symlink: chip-seq/bam → chip-seq/cleanbams"

log "ChIP 3/5: Calling peaks (MACS3)..."
bash "$CHIP_SCRIPTS/03_callpeaks.sh" \
    "$CHIP_DIR" \
    "$REPO_DIR/chip-pipeline/ChIPseq_input_matching.csv" \
    "$PARALLEL_JOBS"
[[ ! -L "$CHIP_DIR/peaks" ]] \
    && ln -s "$CHIP_DIR/macs3" "$CHIP_DIR/peaks" \
    && log "Created symlink: chip-seq/peaks → chip-seq/macs3"

log "ChIP 4a/5: DiffBind spike-in scaling (R)..."
cd "$REPO_DIR" && Rscript "$CHIP_SCRIPTS/04a_diffbind_scale.R" && cd "$INSTALL_DIR"

log "ChIP 4b/5: Generating spike-in bigwigs (bin: ${CHIP_BIGWIG_BIN}bp)..."
bash "$CHIP_SCRIPTS/04b_bamcoverage.sh" \
    "$CHIP_DIR/cleanbams" \
    "$DATA_DIR/ChIPseq_ASCL1_scalefactors.txt" \
    "$THREADS" "$CHIP_BIGWIG_BIN"

log "ChIP 5/5: Peak annotation (R)..."
cd "$REPO_DIR" && Rscript "$CHIP_SCRIPTS/05_peaks_annotation.R" && cd "$INSTALL_DIR"

log "ChIP pipeline complete."
fi


# =============================================================================
# RNA-SEQ PIPELINE
# =============================================================================
if [[ "$RUN_RNA" == true ]]; then

log "══════════════════════════ RNA-SEQ ═══════════════════════════"
eval "$(conda shell.bash hook)" && conda activate atacseq   # kallisto is in atacseq env

log "RNA 1/3: ascl1 experiment (exp1+exp2) — FastQC, TrimGalore, Kallisto..."
bash "$SCRIPT_DIR/run_rnaseq_ascl1.sh"

log "RNA 2/3: hdaci experiment (exp3) — FastQC, TrimGalore, Kallisto..."
bash "$SCRIPT_DIR/run_rnaseq_hdaci.sh"

log "RNA 3/3: phox2b experiment (single-end) — FastQC, TrimGalore, Kallisto..."
bash "$SCRIPT_DIR/run_rnaseq_phox2b.sh"

log "RNA pipeline complete."
fi


log "══════════════════════════════════════════════════"
log "All pipelines complete."
log "Next: render figs/figures.Rmd in RStudio (ascl1_competence/ project root)"
log "══════════════════════════════════════════════════"

RUNSCRIPT

chmod +x "$RUN_SCRIPT"
echo "run_pipeline.sh written to: $RUN_SCRIPT"


# =============================================================================
# DONE
# =============================================================================

echo ""
echo "============================================================"
echo " Setup complete!"
echo "============================================================"
echo ""
echo " Key files:"
echo "   Config    : $CONFIG_FILE"
echo "   Runner    : $RUN_SCRIPT"
echo ""
echo " Directory layout:"
echo "   $INSTALL_DIR/"
echo "   ├── ascl1_competence/   ← repo (config.yaml is inside)"
echo "   ├── Data → $DATA_DIR"
echo "   │   ├── atac-seq/       ← FASTQs + pipeline outputs"
echo "   │   ├── chip-seq/       ← FASTQs + pipeline outputs"
echo "   │   └── rna-seq/        ← FASTQs + pipeline outputs"
echo "   ├── Figures/"
echo "   ├── run_pipeline.sh"
echo "   ├── run_rnaseq_ascl1.sh    ← exp1+exp2 (paired-end)"
echo "   ├── run_rnaseq_hdaci.sh    ← exp3 (paired-end, dynamic sample list)"
echo "   └── run_rnaseq_phox2b.sh  ← phox2b (single-end)"
echo ""
echo " Before running:"
echo "   1. Commit config.yaml to git"
echo ""
echo " To run:"
echo "   bash $RUN_SCRIPT [--atac-only | --chip-only | --rna-only]"
echo "============================================================"
