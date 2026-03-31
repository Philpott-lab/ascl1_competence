#!/bin/bash
set -euo pipefail

# Usage: bash 01_prepare_index.sh <main_directory> [threads]
# Downloads and builds the mm39 bowtie2 index.
# Only needs to be run once per genome.
# threads defaults to 8 if not specified.

if [[ $# -lt 1 ]]; then
    echo "Usage: bash 01_prepare_index.sh <main_directory> [threads]"
    exit 1
fi

MAINDIR=$1
THREADS=${2:-8}
INDEXDIR=$MAINDIR/indexes
mkdir -p $INDEXDIR

echo "Main directory: $MAINDIR"
echo "Threads: $THREADS"

echo "Preparing bowtie2 index..."

# Download genome FASTA
if [[ ! -f "$INDEXDIR/mm39.fa.gz" && ! -f "$INDEXDIR/mm39.fa" ]]; then
    echo "Downloading mm39 genome..."
    wget -P "$INDEXDIR/" http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
else
    echo "Genome FASTA already downloaded. Skipping..."
fi

# Decompress
if [[ ! -f "$INDEXDIR/mm39.fa" ]]; then
    echo "Decompressing genome..."
    gunzip "$INDEXDIR/mm39.fa.gz"
else
    echo "Genome already decompressed. Skipping..."
fi

# Build index
if [[ ! -f "$INDEXDIR/mm39.1.bt2" ]]; then
    echo "Building bowtie2 index..."
    bowtie2-build "$INDEXDIR/mm39.fa" "$INDEXDIR/mm39" --threads "$THREADS"
    echo "Index prepared"
else
    echo "Bowtie2 index already exists. Skipping..."
fi
