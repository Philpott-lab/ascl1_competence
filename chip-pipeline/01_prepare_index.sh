#!/bin/bash
set -euo pipefail

#######################################################################################
##                     ChIP-seq GENOME INDEX PREPARATION (mm39 + dm6)               ##
#######################################################################################


#### OVERVIEW ####
# Downloads mm39 and dm6 genome FASTAs, renames dm6 chromosomes to avoid conflicts,
# concatenates the two genomes, and builds a combined bowtie2 index.
# The dm6 spike-in genome is used for calibration of ChIP-seq signal.
# Only needs to be run once.


#### USAGE ####
# bash 01_prepare_index.sh <main_directory> [threads]
#
# Arguments:
#   $1  MAINDIR  - path to the main project directory
#   $2  THREADS  - threads for bowtie2-build (default: 8)
#
# Example:
#   bash 01_prepare_index.sh /path/to/project 40


#### ARGUMENT HANDLING ####
if [[ $# -lt 1 ]]; then
    echo "Usage: bash 01_prepare_index.sh <main_directory> [threads]"
    exit 1
fi

MAINDIR=$1
THREADS=${2:-8}

INDEXDIR=$MAINDIR/indexes
mkdir -p "$INDEXDIR"

MM39_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz"
DM6_URL="http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz"

MM39_FA="$INDEXDIR/mm39.fa"
DM6_FA="$INDEXDIR/dm6.fa"
RENAMED_DM6_FA="$INDEXDIR/dm6_renameID.fa"
COMBINED_FA="$INDEXDIR/mm39_dm6.fa"
INDEX_PREFIX="$INDEXDIR/mm39_dm6"

echo "Index directory : $INDEXDIR"
echo "Threads         : $THREADS"


#### CHECK IF INDEX ALREADY EXISTS ####
if [[ -f "${INDEX_PREFIX}.1.bt2" ]]; then
    echo "Bowtie2 index already exists at ${INDEX_PREFIX}. Skipping."
    exit 0
fi


#######################################################################################
##                              Download genome FASTAs                               ##
#######################################################################################

# Download mm39
if [[ ! -f "$MM39_FA" && ! -f "${MM39_FA}.gz" ]]; then
    echo "Downloading mm39 genome..."
    wget -nc -P "$INDEXDIR/" "$MM39_URL"
else
    echo "mm39 FASTA already downloaded. Skipping..."
fi

# Download dm6
if [[ ! -f "$DM6_FA" && ! -f "${DM6_FA}.gz" ]]; then
    echo "Downloading dm6 genome..."
    wget -nc -P "$INDEXDIR/" "$DM6_URL"
else
    echo "dm6 FASTA already downloaded. Skipping..."
fi


#######################################################################################
##                              Decompress FASTAs                                    ##
#######################################################################################

if [[ ! -f "$MM39_FA" ]]; then
    echo "Decompressing mm39..."
    gunzip "${MM39_FA}.gz"
else
    echo "mm39 FASTA already decompressed. Skipping..."
fi

if [[ ! -f "$DM6_FA" ]]; then
    echo "Decompressing dm6..."
    gunzip "${DM6_FA}.gz"
else
    echo "dm6 FASTA already decompressed. Skipping..."
fi


#######################################################################################
##                    Rename dm6 chromosomes and concatenate                         ##
#######################################################################################

# Rename dm6 chromosomes from e.g. chrX to dm6_chrX to avoid name collisions with mm39
# This allows reads to be separated cleanly by chromosome name prefix after alignment
if [[ ! -f "$RENAMED_DM6_FA" ]]; then
    echo "Renaming dm6 chromosomes (chr -> dm6_chr)..."
    sed 's/chr/dm6_chr/g' "$DM6_FA" > "$RENAMED_DM6_FA"
else
    echo "Renamed dm6 FASTA already exists. Skipping..."
fi

if [[ ! -f "$COMBINED_FA" ]]; then
    echo "Concatenating mm39 and dm6 genomes..."
    cat "$MM39_FA" "$RENAMED_DM6_FA" > "$COMBINED_FA"
else
    echo "Combined FASTA already exists. Skipping..."
fi


#######################################################################################
##                              Build bowtie2 index                                  ##
#######################################################################################

echo "Building bowtie2 index (this may take a while)..."
bowtie2-build "$COMBINED_FA" "$INDEX_PREFIX" --threads "$THREADS" --quiet

echo "Index prepared at: $INDEX_PREFIX"
