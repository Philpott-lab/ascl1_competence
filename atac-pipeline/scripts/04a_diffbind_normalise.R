#####################################################################################
##           ATAC-seq PEAK COUNTING AND TMM NORMALISATION WITH DIFFBIND           ##
#####################################################################################


#### OVERVIEW ####
# 1. Loads sample metadata and MACS3 peak files into a DiffBind object
# 2. Applies mm39 blacklist regions (ENCODE + excluderanges)
# 3. Builds a consensus peakset: peaks present in >= 3/4 samples per condition
# 4. Filters to canonical chromosomes (chr1-19, X, Y)
# 5. Counts reads across the consensus peakset
# 6. Applies TMM normalisation and exports scale factors for bamCoverage


#### INPUT DATA ####
# Samples from mESCs with Ascl1-GR-HA
# 3 conditions: mESC, EpiLC, NE
# 3 timepoints per condition: 0h, 6h, 24h after dexamethasone induction of Ascl1
# DiffBind sample sheet expected at: ../Data/atac-seq/atacseq_diffbindSamples.csv
# Peak files expected as produced by 03_call_peaks.sh (MACS3, narrowPeak format)


#### PACKAGES ####
# Install via Bioconductor:
#   BiocManager::install(c("DiffBind", "AnnotationHub", "GenomicRanges"))
# Install via CRAN:
#   install.packages(c("here", "tidyverse", "ggplot2", "cowplot", "edgeR"))


#### EXECUTING THE SCRIPT ####
# Run interactively in RStudio or from the command line:
#   Rscript 04a_diffbind_count_normalise.R
# Outputs:
#   ../Data/ATACseq_counted_nochrUn.Rds   — counted DiffBind object
#   ../Data/ATACseq_analysed_TMM.Rds      — TMM-normalised DiffBind object
#   ../Data/ATAC_TMM_scalefactors.csv     — TMM scale factors for bamCoverage


#####################################################################################

library(here)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(DiffBind)
library(edgeR)
library(AnnotationHub)
library(GenomicRanges)


#####################################################################################
##                    1. Build DiffBind object and apply blacklist                 ##
#####################################################################################

atac_samples <- read.csv(here("../Data/atac-seq/atacseq_diffbindSamples.csv"))

# mm39 blacklist: DiffBind's built-in blacklist uses ENCODE regions but mm39 is not
# yet available there. Instead we use the excluderanges resource via AnnotationHub.
# mm39 blacklist ID: AH107321
# Source: https://github.com/dozmorovlab/excluderanges
ah <- AnnotationHub()
atac_blacklist_mm39 <- ah[["AH107321"]]
rm(ah)

# Build initial DiffBind object (minOverlap = 1 retains all peaks at this stage)
atac_data <- dba(sampleSheet = atac_samples, minOverlap = 1)

# Apply blacklist — greylist = FALSE as we have no input/control tracks
atac_data <- dba.blacklist(atac_data,
                           blacklist = atac_blacklist_mm39,
                           greylist  = FALSE)
# Expected output:
#   Removed: 198227 of 4703301 intervals.
#   Removed: 13827 merged (of 355316) and 81704 (of 355316) consensus.


#####################################################################################
##                       2. Build consensus peakset                                ##
#####################################################################################

# Consensus per condition: keep peaks present in >= 3 out of 4 replicates
# This reduces noise from low-confidence peaks while retaining condition-specific sites
atac_consensus_perCond <- dba.peakset(atac_data,
                                      consensus   = DBA_TREATMENT,
                                      minOverlap  = 3)

# Merge into a single consensus peakset across all conditions
atac_consensus_full <- dba(atac_consensus_perCond,
                           mask       = atac_consensus_perCond$masks$Consensus,
                           minOverlap = 1)

# Retrieve as a data frame for filtering
atac_consensus_peakset <- dba.peakset(atac_consensus_full,
                                      bRetrieve  = TRUE,
                                      DataType   = DBA_DATA_FRAME,
                                      minOverlap = 1)

# Filter to canonical chromosomes only (chr1-19, X, Y)
# Removes unplaced/unlocalized contigs (chrUn_*, chr*_random, etc.)
# Expected: removes ~93 peaks (202,644 -> 202,551)
atac_consensus_peakset_filt <- atac_consensus_peakset[
  grepl("chr([1-9]|1[0-9]|X|Y)$", atac_consensus_peakset$CHR), ]

message("Chromosomes retained: ", paste(unique(atac_consensus_peakset_filt$CHR), collapse = ", "))
message("Peaks after filtering: ", nrow(atac_consensus_peakset_filt))


#####################################################################################
##                         3. Count reads across consensus peaks                   ##
#####################################################################################

# summits = 75: re-centres each peak and extends 75bp either side (150bp windows)
# This value follows Rory Stark's recommendation that most ATAC fragments are 50-100bp
# bParallel = TRUE: uses BiocParallel for faster counting
atac_dba_counted <- dba.count(atac_data,
                              peaks     = atac_consensus_peakset_filt,
                              summits   = 75,
                              bParallel = TRUE)

saveRDS(atac_dba_counted, file = here("../Data/ATACseq_counted_nochrUn.Rds"))
message("Counted DiffBind object saved.")


#####################################################################################
##                      4. TMM normalisation and scale factor export               ##
#####################################################################################

# Load counted object (allows this section to be re-run independently)
atac_dba <- readRDS(file = here("../Data/ATACseq_counted_nochrUn.Rds"))

# Normalise using TMM on peak reads (RiP), no background correction
# DBA_LIBSIZE_PEAKREADS: library size = reads in peaks only (appropriate for ATAC)
# background = FALSE: no background bin normalisation
atac_normalised <- dba.normalize(atac_dba,
                                 normalize = DBA_NORM_TMM,
                                 library   = DBA_LIBSIZE_PEAKREADS,
                                 background = FALSE)

saveRDS(atac_normalised, file = here("../Data/ATACseq_analysed_TMM.Rds"))
message("Normalised DiffBind object saved.")

# Retrieve normalisation factors
norm_info <- dba.normalize(atac_normalised, bRetrieve = TRUE)

# Build scale factor table
# bam_coverage_scale = 1 / tmm_scale_factor
# bamCoverage multiplies coverage by --scaleFactor, so the TMM factor must be inverted:
# a library that is 2x larger gets a TMM factor of ~2 and a bamCoverage scale of ~0.5
scale_factor_df <- data.frame(
  sample             = names(norm_info$norm.factors),
  lib_size           = norm_info$lib.sizes,
  tmm_scale_factor   = norm_info$norm.factors,
  bam_coverage_scale = 1 / norm_info$norm.factors,
  row.names          = NULL
)

print(scale_factor_df)

write.csv(scale_factor_df,
          file      = here("../Data/ATAC_TMM_scalefactors.csv"),
          row.names = FALSE)

message("Scale factors written to ../Data/ATAC_TMM_scalefactors.csv")
message("Pass bam_coverage_scale values to bamCoverage --scaleFactor in 04b_bigwig_tmm.sh")
