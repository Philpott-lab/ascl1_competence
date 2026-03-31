#####################################################################################
##          ChIP-seq PEAK SET GENERATION, BED FILE EXPORT AND ANNOTATION          ##
#####################################################################################


#### OVERVIEW ####
# 1. Extracts top 10,000 peaks per condition ranked by average read intensity
# 2. Writes per-condition BED files and HOMER-formatted BED files
# 3. Computes overlaps between condition peak sets and assigns peaks to categories:
#    universal (all 3), condition-specific, and pairwise shared
# 4. Writes overlap category BED files and HOMER-formatted versions
# 5. Finds genes within 200kb of each peak set (per condition and per overlap category)
# 6. Writes peak-to-gene association tables


#### INPUT DATA ####
# Normalised DiffBind object: ../Data/ChIPseq_ASCL1_normalised.Rds
# (produced by chip_05a_diffbind_count_normalise.R)


#### PACKAGES ####
# Install via Bioconductor:
#   BiocManager::install(c("DiffBind", "GenomicRanges", "rtracklayer",
#                          "AnnotationHub", "BSgenome.Mmusculus.UCSC.mm39",
#                          "GenomicInteractions"))
# Install via CRAN:
#   install.packages(c("here", "tidyverse", "eulerr"))


#### EXECUTING THE SCRIPT ####
# Run interactively in RStudio or from the command line:
#   Rscript chip_06_peaksets_annotation.R
# Outputs (all in ../Data/):
#   ChIPseq_ASCL1_mesc10k.bed / _epilc10k.bed / _ne10k.bed
#   ChIPseq_ASCL1_mesc10k_homer.bed / _epilc10k_homer.bed / _ne10k_homer.bed
#   ChIPseq_ASCL1_universal10k.bed
#   ChIPseq_ASCL1_mesconly10k.bed / _epilconly10k.bed / _neonly10k.bed
#   ChIPseq_ASCL1_mescepilc10k.bed / _mescne10k.bed / _epilcne10k.bed
#   ChIPseq_ASCL1_universal10k_homer.bed
#   ChIPseq_ASCL1_mesconly10k_homer.bed / _epilconly10k_homer.bed / _neonly10k_homer.bed
#   ChIP_ASCL1_mESC_peakstogenes.txt / _EpiLC_ / _NE_
#   ChIPseq_allpeaktogenes.txt


#####################################################################################

library(here)
library(tidyverse)
library(DiffBind)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm39)
library(GenomicInteractions)
library(eulerr)


#####################################################################################
##                     1. Extract top 10,000 peaks per condition                   ##
#####################################################################################

chip_dba     <- readRDS(here("../Data/ChIPseq_ASCL1_normalised.Rds"))
chip_peakset <- dba.peakset(chip_dba, bRetrieve = TRUE)

# Calculate average read intensity across replicates per condition
chip_peakset$avg_gG5_mESC  <- rowMeans(as.data.frame(mcols(chip_peakset)[, grep("gG5_mESC",  colnames(mcols(chip_peakset)))]))
chip_peakset$avg_gG5_EpiLC <- rowMeans(as.data.frame(mcols(chip_peakset)[, grep("gG5_EpiLC", colnames(mcols(chip_peakset)))]))
chip_peakset$avg_gG5_NE    <- rowMeans(as.data.frame(mcols(chip_peakset)[, grep("gG5_NE",    colnames(mcols(chip_peakset)))]))

# Extract top 10,000 sites per condition ranked by average intensity
top10k_mESC  <- chip_peakset[order(chip_peakset$avg_gG5_mESC,  decreasing = TRUE)][1:10000] %>% as.data.frame()
top10k_EpiLC <- chip_peakset[order(chip_peakset$avg_gG5_EpiLC, decreasing = TRUE)][1:10000] %>% as.data.frame()
top10k_NE    <- chip_peakset[order(chip_peakset$avg_gG5_NE,    decreasing = TRUE)][1:10000] %>% as.data.frame()

message("Top 10k peaks extracted per condition")


#####################################################################################
##                    2. Write per-condition BED and HOMER BED files               ##
#####################################################################################

# --- Standard BED files (chr, start, end) ---

write.table(top10k_mESC[,  c(1,2,3)],
            file = here("../Data/ChIPseq_ASCL1_mesc10k.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(top10k_EpiLC[, c(1,2,3)],
            file = here("../Data/ChIPseq_ASCL1_epilc10k.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(top10k_NE[,    c(1,2,3)],
            file = here("../Data/ChIPseq_ASCL1_ne10k.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# --- HOMER-formatted BED files (chr, start, end, ID, NA, strand) ---
# HOMER requires 6 columns: chr, start, end, unique peak ID, score (NA), strand

write_homer_bed <- function(bed_df, prefix, outfile) {
  read_tsv(bed_df, col_names = c("chr", "start", "end"),
           show_col_types = FALSE) %>%
    mutate(ID     = paste0(prefix, "_peak_", row_number()),
           na     = "NA",
           strand = "*") %>%
    write_tsv(outfile, col_names = FALSE)
}

write_homer_bed(here("../Data/ChIPseq_ASCL1_mesc10k.bed"),
                "mesc",  here("../Data/ChIPseq_ASCL1_mesc10k_homer.bed"))
write_homer_bed(here("../Data/ChIPseq_ASCL1_epilc10k.bed"),
                "epilc", here("../Data/ChIPseq_ASCL1_epilc10k_homer.bed"))
write_homer_bed(here("../Data/ChIPseq_ASCL1_ne10k.bed"),
                "ne",    here("../Data/ChIPseq_ASCL1_ne10k_homer.bed"))

message("Per-condition BED and HOMER BED files written")


#####################################################################################
##              3. Compute overlaps and assign peaks to categories                 ##
#####################################################################################

# Load per-condition peak sets as GRanges
mesc  <- read_table(here("../Data/ChIPseq_ASCL1_mesc10k.bed"),
                    col_names = FALSE, show_col_types = FALSE) %>%
  setNames(c("chr", "start", "end", "peakID")) %>% GRanges()
epilc <- read_table(here("../Data/ChIPseq_ASCL1_epilc10k.bed"),
                    col_names = FALSE, show_col_types = FALSE) %>%
  setNames(c("chr", "start", "end", "peakID")) %>% GRanges()
ne    <- read_table(here("../Data/ChIPseq_ASCL1_ne10k.bed"),
                    col_names = FALSE, show_col_types = FALSE) %>%
  setNames(c("chr", "start", "end", "peakID")) %>% GRanges()

names(mesc)  <- paste0("mesc_",  seq_along(mesc))
names(epilc) <- paste0("epilc_", seq_along(epilc))
names(ne)    <- paste0("ne_",    seq_along(ne))

# Combine all peaks and reduce to non-overlapping intervals
all_peaks        <- c(mesc, epilc, ne)
all_peaks$source <- c(rep("mesc", length(mesc)),
                      rep("epilc", length(epilc)),
                      rep("ne", length(ne)))

reduced_peaks <- reduce(all_peaks, with.revmap = TRUE)

# Label each reduced peak by which conditions it originated from
peak_membership <- lapply(reduced_peaks$revmap, function(ix) {
  sources <- all_peaks[ix]$source
  paste0(sort(unique(sources)), collapse = "_")
})
reduced_peaks$group <- unlist(peak_membership)

# Extract overlap categories
universal  <- reduced_peaks[reduced_peaks$group == "epilc_mesc_ne"]
mesc_only  <- reduced_peaks[reduced_peaks$group == "mesc"]
epilc_only <- reduced_peaks[reduced_peaks$group == "epilc"]
ne_only    <- reduced_peaks[reduced_peaks$group == "ne"]
mesc_epilc <- reduced_peaks[reduced_peaks$group == "epilc_mesc"]
mesc_ne    <- reduced_peaks[reduced_peaks$group == "mesc_ne"]
epilc_ne   <- reduced_peaks[reduced_peaks$group == "epilc_ne"]

message("Overlap categories:")
message("  Universal (all 3):  ", length(universal))
message("  mESC only:          ", length(mesc_only))
message("  EpiLC only:         ", length(epilc_only))
message("  NE only:            ", length(ne_only))
message("  mESC + EpiLC:       ", length(mesc_epilc))
message("  mESC + NE:          ", length(mesc_ne))
message("  EpiLC + NE:         ", length(epilc_ne))


#####################################################################################
##                    4. Write overlap category BED and HOMER BED files            ##
#####################################################################################

write_bed <- function(gr, outfile) {
  write.table(as.data.frame(gr)[, c(1,2,3)],
              file = outfile,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

write_bed(universal,  here("../Data/ChIPseq_ASCL1_universal10k.bed"))
write_bed(mesc_only,  here("../Data/ChIPseq_ASCL1_mesconly10k.bed"))
write_bed(epilc_only, here("../Data/ChIPseq_ASCL1_epilconly10k.bed"))
write_bed(ne_only,    here("../Data/ChIPseq_ASCL1_neonly10k.bed"))
write_bed(mesc_epilc, here("../Data/ChIPseq_ASCL1_mescepilc10k.bed"))
write_bed(mesc_ne,    here("../Data/ChIPseq_ASCL1_mescne10k.bed"))
write_bed(epilc_ne,   here("../Data/ChIPseq_ASCL1_epilcne10k.bed"))

# HOMER-formatted versions of the overlap category BED files
write_homer_bed(here("../Data/ChIPseq_ASCL1_universal10k.bed"),
                "mesc",  here("../Data/ChIPseq_ASCL1_universal10k_homer.bed"))
write_homer_bed(here("../Data/ChIPseq_ASCL1_mesconly10k.bed"),
                "mesc",  here("../Data/ChIPseq_ASCL1_mesconly10k_homer.bed"))
write_homer_bed(here("../Data/ChIPseq_ASCL1_epilconly10k.bed"),
                "epilc", here("../Data/ChIPseq_ASCL1_epilconly10k_homer.bed"))
write_homer_bed(here("../Data/ChIPseq_ASCL1_neonly10k.bed"),
                "ne",    here("../Data/ChIPseq_ASCL1_neonly10k_homer.bed"))

message("Overlap category BED and HOMER BED files written")


#####################################################################################
##              5. Find genes within 200kb of per-condition peak sets              ##
#####################################################################################

# Load Ensembl mm39 gene annotation via AnnotationHub
# EnsDb.Mmusculus.v108 AnnotationHub ID: AH109367
ah <- AnnotationHub()
EnsDb.Mmusculus.v108 <- ah[["AH109367"]]
ensembldb::seqlevelsStyle(EnsDb.Mmusculus.v108) <- "UCSC"
genes_gr <- genes(EnsDb.Mmusculus.v108, return.type = "GRanges")
rm(ah)

# Helper function: find all genes within a window around each peak
find_nearby_genes <- function(peaks_gr, genes_gr, window = 200000, set_name = "") {
  peak_results <- lapply(seq_along(peaks_gr), function(i) {
    peak        <- peaks_gr[i]
    peak_window <- resize(peak, width = window, fix = "center")
    overlaps    <- findOverlaps(genes_gr, peak_window)
    if (length(overlaps) == 0) return(NULL)
    nearby_genes <- genes_gr[queryHits(overlaps)]
    data.frame(
      peak_index  = i,
      chr         = as.character(seqnames(peak)),
      peak_start  = start(peak),
      peak_end    = end(peak),
      gene_id     = mcols(nearby_genes)$gene_id,
      gene_symbol = mcols(nearby_genes)$symbol,
      peakset     = set_name,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, peak_results)
}

# Load per-condition top 10k peak sets as GRanges
top10k_gr <- list(
  mESC  = read.table(here("../Data/ChIPseq_ASCL1_mesc10k.bed"))  %>%
    dplyr::rename(chr = V1, start = V2, end = V3) %>% GRanges(),
  EpiLC = read.table(here("../Data/ChIPseq_ASCL1_epilc10k.bed")) %>%
    dplyr::rename(chr = V1, start = V2, end = V3) %>% GRanges(),
  NE    = read.table(here("../Data/ChIPseq_ASCL1_ne10k.bed"))    %>%
    dplyr::rename(chr = V1, start = V2, end = V3) %>% GRanges()
)

# Run per-condition peak-to-gene association
results_per_condition <- lapply(names(top10k_gr), function(set_name) {
  message("Finding nearby genes for: ", set_name)
  find_nearby_genes(top10k_gr[[set_name]], genes_gr, set_name = set_name)
})
names(results_per_condition) <- names(top10k_gr)

write.table(results_per_condition$mESC,
            file = here("../Data/ChIP_ASCL1_mESC_peakstogenes.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(results_per_condition$EpiLC,
            file = here("../Data/ChIP_ASCL1_EpiLC_peakstogenes.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(results_per_condition$NE,
            file = here("../Data/ChIP_ASCL1_NE_peakstogenes.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

message("Per-condition peak-to-gene tables written")


#####################################################################################
##           6. Find genes within 200kb of overlap category peak sets              ##
#####################################################################################

overlap_gr <- list(
  mESC      = read_delim(here("../Data/ChIPseq_ASCL1_mesconly10k.bed"),
                         col_names = FALSE, show_col_types = FALSE) %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    makeGRangesFromDataFrame(keep.extra.columns = FALSE),
  mESC_EpiLC = read_delim(here("../Data/ChIPseq_ASCL1_mescepilc10k.bed"),
                           col_names = FALSE, show_col_types = FALSE) %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    makeGRangesFromDataFrame(keep.extra.columns = FALSE),
  EpiLC     = read_delim(here("../Data/ChIPseq_ASCL1_epilconly10k.bed"),
                         col_names = FALSE, show_col_types = FALSE) %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    makeGRangesFromDataFrame(keep.extra.columns = FALSE),
  EpiLC_NE  = read_delim(here("../Data/ChIPseq_ASCL1_epilcne10k.bed"),
                          col_names = FALSE, show_col_types = FALSE) %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    makeGRangesFromDataFrame(keep.extra.columns = FALSE),
  NE        = read_delim(here("../Data/ChIPseq_ASCL1_neonly10k.bed"),
                         col_names = FALSE, show_col_types = FALSE) %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    makeGRangesFromDataFrame(keep.extra.columns = FALSE),
  mESC_NE   = read_delim(here("../Data/ChIPseq_ASCL1_mescne10k.bed"),
                          col_names = FALSE, show_col_types = FALSE) %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    makeGRangesFromDataFrame(keep.extra.columns = FALSE),
  Universal = read_delim(here("../Data/ChIPseq_ASCL1_universal10k.bed"),
                          col_names = FALSE, show_col_types = FALSE) %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    makeGRangesFromDataFrame(keep.extra.columns = FALSE)
)

results_overlap <- lapply(names(overlap_gr), function(set_name) {
  message("Finding nearby genes for overlap category: ", set_name)
  find_nearby_genes(overlap_gr[[set_name]], genes_gr, set_name = set_name)
})
names(results_overlap) <- names(overlap_gr)

# Combine all overlap categories into a single table
final_results <- do.call(rbind, results_overlap)
message("Total peak-gene interactions across overlap categories: ", nrow(final_results))

write.table(final_results,
            file = here("../Data/ChIPseq_allpeaktogenes.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

message("All peak-to-gene annotation tables written.")
