Input files required
RNAseq (Experiment 1)

../Data/RNAseq_ASCL1_dds.Rds — DESeq2 object
../Data/RNAseq_ASCL1_DEGs.csv — DEG table with direction/timepoint columns
../Data/RNAseq_ASCL1_zmtrx.txt — z-scored expression matrix

ChIP-seq (Experiment 1)

../Data/ChIPseq_ASCL1_samples.csv — DiffBind sample sheet
../Data/ChIPseq_ASCL1_normalised.Rds — normalised DiffBind object
../Data/ChIPseq_ASCL1_{mesc/epilc/ne}10k.bed — consensus peaksets per condition
../Data/ChIPseq_ASCL1_{mesconly/mescepilc/epilconly/epilcne/neonly/mescne/universal}10k.bed — split peaksets
../Data/ChIPseq_allpeaktogenes.txt — peak-to-gene mapping table
../Data/ChIPseq_NEspecificDEG_peaks_full.csv — peaks near NE-specific DEGs
../Data/chipseq_ascl1-ha/bigwig/*.bigwig — merged bigWig files
../Data/chipseq_published/ — published ChIP-seq BED files (borromeo2014, webb2013, raposo2015, casey2018, wapinski2013, wang2023)
../Data/chipseq_published/cruzmolina2017/ — H3K27ac bigWigs (2i, AntNPC)
../Data/chipseq_published/bleckwehl2020/ — H3K4me1 bigWig (mESC)
../Data/homer/ne_specific/knownResults.txt — HOMER motif results (NE)
../Data/homer/mesc_specific/knownResults.txt — HOMER motif results (mESC)

ATAC-seq

../Data/ATACseq_analysed_RLE.Rds — analysed DiffBind object
../Data/ATACseq/atacseq_diffbindSamples.csv — DiffBind sample sheet
../Data/ATACseq_counted_nochrUn.Rds — counted DiffBind object (for FRiP)
../Data/ATACseq_{mESC/EpiLC/NE}_{0h/6h/24h}_peaks_in3of4.bed — consensus ATAC peaks
../Data/ATACseq/tss_enrichment/readCountInPeaks.txt.summary
../Data/ATACseq/tss_enrichment/readCountInFlanks.txt.summary
../Data/ATACseq/mtDNA_content.txt
../Data/atacseq/bigwig/*.bigwig — merged ATAC bigWig files

RNAseq (Experiments 2 & 3 — HDACi / Cofactor)

../Data/RNAseq_HDACi_dds.rds
../Data/RNAseq_HDACi_normcounts.csv
../Data/RNAseq_HDACi_DEGs.csv
../Data/RNAseq_HDACi_DEGs_mESC.csv
../Data/RNAseq_GeneMap.csv
../Data/RNAseq_Cofactor_Phox2b_DEGs.csv
../Data/RNAseq_Cofactors_Phox2b_normcounts.csv

Reference / annotation

../Data/Ref_MouseTFs.csv
Nuclei_2.csv

------

Intermediate files produced by the figures script: 

../Data/ChIPseq_NEboundgenes_208.txt
../Data/ChIPseq_653NEsites.bed + split versions:

ChIPseq_653NEsites_splitNEonly.bed
ChIPseq_653NEsites_splitnotNEonly.bed
ChIPseq_653NEsites_splitnotmesc.bed
ChIPseq_653NEsites_splitmesc.bed
ChIPseq_653NEsites_closedT0.bed
ChIPseq_653NEsites_openT0.bed
ChIPseq_653NEsites_opened.bed


../Data/ChIPseq_362_openmesc_boundmesc.bed
../Data/ChIPseq_362_openmesc_notboundmesc.bed
../Data/cluster1_genes.csv, cluster2_genes.csv
../Data/cluster1_GO.csv, cluster2_GO.csv
../Data/cluster_phox2b_lower_genes.csv, cluster_phox2b_lower_GO.csv
deeptools .mat.gz matrix files (many, spread across bash chunks)

