#####################################################################################
##        ChIP-seq PEAK COUNTING, NORMALISATION AND SCALE FACTOR EXPORT           ##
#####################################################################################


#### OVERVIEW ####
# 1. Builds a custom blacklist combining mm39 ENCODE regions and C11 control peaks
# 2. Loads sample metadata and MACS3 peak files into a DiffBind object
# 3. Builds a consensus peakset: peaks present in >= 2/4 replicates per condition
# 4. Counts reads across the consensus peakset
# 5. Applies RLE normalisation on peak reads and runs differential analysis
# 6. Exports RLE scale factors for use in bamCoverage (script 05b)


#### INPUT DATA ####
# Samples from HA-Ascl1 ChIP in three conditions: mESC, EpiLC, NE
# C11 cells used as negative control (no Ascl1 expression) — peaks used as blacklist
# DiffBind sample sheet expected at: ../Data/config/chip_DiffBind_samples.csv
#   Tissue column should contain "gG5" for experimental samples and "C11" for controls
# Peak files expected as produced by chip_03_call_peaks.sh (MACS3, narrowPeak format)
# MACS3 peak directory expected at: ../Data/macs3/


#### PACKAGES ####
# Install via Bioconductor:
#   BiocManager::install(c("DiffBind", "AnnotationHub", "GenomicRanges",
#                          "BSgenome.Mmusculus.UCSC.mm39"))
# Install via CRAN:
#   install.packages(c("here", "tidyverse"))


#### EXECUTING THE SCRIPT ####
# Run interactively in RStudio or from the command line:
#   Rscript chip_05a_diffbind_count_normalise.R
# Outputs:
#   ../Data/ChIPseq_counted_ascl1_peaks.Rds      — counted DiffBind object
#   ../Data/ChIPseq_ASCL1_normalised.Rds         — RLE-normalised DiffBind object
#   ../Data/ChIPseq_ASCL1_scalefactors.txt        — RLE scale factors for bamCoverage


#####################################################################################

library(here)
library(tidyverse)
library(DiffBind)
library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm39)


#####################################################################################
##           1. Build custom blacklist: ENCODE mm39 + C11 control peaks            ##
#####################################################################################

genome_mm39 <- BSgenome.Mmusculus.UCSC.mm39

# Load C11 control peaks as an additional blacklist
# C11 cells do not express Ascl1, so any peaks called in this condition represent
# non-specific signal and should be excluded from downstream analysis
c11_peak_files <- list.files(
  path    = here("../Data/macs3/"),
  pattern = "^C11.*\\.narrowPeak$",
  full.names = FALSE
)

if (length(c11_peak_files) == 0) {
  stop("[ERROR]: No C11 narrowPeak files found in ../Data/macs3/. Check peak directory.")
}

c11_list <- lapply(c11_peak_files, function(file) {
  read.table(file.path(here("../Data/macs3/"), file), header = FALSE, sep = "\t")
})

c11_combined    <- do.call(rbind, c11_list)
c11_combined_gr <- GRanges(
  seqnames = c11_combined$V1,
  ranges   = IRanges(start = c11_combined$V2, end = c11_combined$V3),
  name     = "C11 Control"
)
seqlengths(c11_combined_gr) <- seqlengths(genome_mm39)[seqlevels(c11_combined_gr)]
genome(c11_combined_gr)     <- "mm39"

message("C11 control peaks loaded: ", length(c11_combined_gr), " regions")

# Load mm39 ENCODE blacklist via AnnotationHub
# mm39 is not yet available in DiffBind's built-in blacklist so we use excluderanges
# Source: https://github.com/dozmorovlab/excluderanges
# mm39 AnnotationHub ID: AH107321
ah <- AnnotationHub()
chip_blacklist_mm39        <- ah[["AH107321"]]
genome(chip_blacklist_mm39) <- "mm39"
rm(ah)

# Remove any C11 peaks that overlap the ENCODE blacklist (already covered)
# then combine with ENCODE blacklist for a single unified exclusion set
overlaps              <- findOverlaps(c11_combined_gr, chip_blacklist_mm39)
c11_filtered          <- c11_combined_gr[-queryHits(overlaps)]
combined_blacklist    <- c(c11_filtered, chip_blacklist_mm39)
combined_blacklist    <- trim(combined_blacklist)

message("Combined blacklist regions: ", length(combined_blacklist))


#####################################################################################
##                    2. Build DiffBind object and apply blacklist                 ##
#####################################################################################

# Set working directory for relative BAM paths in the DiffBind sample sheet
setwd(here(".."))

chip_samples <- read.csv(here("../Data/config/chip_DiffBind_samples.csv"))

# Retain only experimental gG5 samples; C11 samples used only for blacklist above
chip_samples <- chip_samples %>% dplyr::filter(Tissue == "gG5")

# Build initial DiffBind object (minOverlap = 1 retains all peaks at this stage)
chip_dba <- dba(sampleSheet = chip_samples, minOverlap = 1)

# Apply combined blacklist
chip_dba <- dba.blacklist(chip_dba,
                          blacklist = combined_blacklist,
                          greylist  = FALSE)


#####################################################################################
##                         3. Build consensus peakset                              ##
#####################################################################################

# Consensus per condition (Tissue + Condition):
# keep peaks present in >= 2 out of 4 replicates per condition group
chip_consensus_perCond <- dba.peakset(chip_dba,
                                      consensus  = c(DBA_TISSUE, DBA_CONDITION),
                                      minOverlap = 2)

# Merge into a single consensus peakset across all conditions
chip_consensus_full <- dba(chip_consensus_perCond,
                           mask       = chip_consensus_perCond$masks$Consensus,
                           minOverlap = 1)

# Retrieve as data frame for counting
chip_consensus_peakset <- dba.peakset(chip_consensus_full,
                                      bRetrieve = TRUE,
                                      DataType  = DBA_DATA_FRAME)

message("Consensus peaks: ", nrow(chip_consensus_peakset))


#####################################################################################
##                         4. Count reads across consensus peaks                   ##
#####################################################################################

# Check if counted object already exists to avoid re-running the slow counting step
counted_path <- here("../Data/ChIPseq_counted_ascl1_peaks.Rds")

if (!file.exists(counted_path)) {
  message("Counting reads across consensus peakset...")
  # summits = 250: re-centres each peak and extends 250bp either side (500bp windows)
  # bUseSummarizeOverlaps = TRUE: more accurate counting for ChIP-seq
  chip_dba_counted <- dba.count(chip_dba,
                                peaks                = chip_consensus_peakset,
                                summits              = 250,
                                bParallel            = TRUE,
                                bUseSummarizeOverlaps = TRUE)
  saveRDS(chip_dba_counted, file = counted_path)
  message("Counted DiffBind object saved to: ", counted_path)
} else {
  chip_dba_counted <- readRDS(counted_path)
  message("Loaded existing counted DiffBind object from: ", counted_path)
}


#####################################################################################
##                    5. RLE normalisation and differential analysis               ##
#####################################################################################

# Load counted object (allows this section to be re-run independently)
chip_dba <- readRDS(counted_path)

# Normalise using RLE (native DESeq2 method) on peak reads
# DBA_NORM_NATIVE with DESeq2/edgeR uses RLE and TMM respectively
# DBA_LIBSIZE_PEAKREADS: library size = reads in peaks only
chip_dba <- dba.normalize(chip_dba,
                          method    = DBA_ALL_METHODS,
                          normalize = DBA_NORM_NATIVE,
                          library   = DBA_LIBSIZE_PEAKREADS)

chip_dba <- dba.analyze(chip_dba, method = DBA_ALL_METHODS)

saveRDS(chip_dba, file = here("../Data/ChIPseq_ASCL1_normalised.Rds"))
message("Normalised DiffBind object saved.")


#####################################################################################
##                         6. Export RLE scale factors                             ##
#####################################################################################

# Load normalised object (allows this section to be re-run independently)
chip_dba <- readRDS(here("../Data/ChIPseq_ASCL1_normalised.Rds"))

# Retrieve the RLE size factors from the DESeq2 analysis
norm <- dba.analyze(chip_dba, bRetrieveAnalysis = TRUE)

# Build scale factor table
# InvNormFacs = 1 / NormFacs is the value passed to bamCoverage --scaleFactor
# A library with a size factor > 1 (larger than average) gets scaled down (< 1)
scale_factor_df <- data.frame(
  bamID        = c("gG5_mESC_HA_1_markdup.bam",
                   "gG5_mESC_HA_2_markdup.bam",
                   "gG5_mESC_HA_3_markdup.bam",
                   "gG5_mESC_HA_4_markdup.bam",
                   "gG5_EpiLC_HA_1_markdup.bam",
                   "gG5_EpiLC_HA_2_markdup.bam",
                   "gG5_EpiLC_HA_3_markdup.bam",
                   "gG5_EpiLC_HA_4_markdup.bam",
                   "gG5_NE_HA_1_markdup.bam",
                   "gG5_NE_HA_2_markdup.bam",
                   "gG5_NE_HA_4_markdup.bam"),
  NormFacs     = norm$sizeFactor,
  InvNormFacs  = 1 / norm$sizeFactor
)

print(scale_factor_df)

write_delim(scale_factor_df,
            file = here("../Data/ChIPseq_ASCL1_scalefactors.txt"))

message("Scale factors written to ../Data/ChIPseq_ASCL1_scalefactors.txt")
message("Pass InvNormFacs values to bamCoverage --scaleFactor in chip_05b_bigwig.sh")
