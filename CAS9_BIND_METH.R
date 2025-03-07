# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor and CRAN packages
#BiocManager::install(c("Rsubread", "DESeq2", "ggplot2", "tidyverse", "VennDiagram"), ask = FALSE)

# Load libraries
library(VennDiagram)
library(Rsubread)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(data.table)
library(dplyr)

base_dir<- "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore"

cas9_genes_all <- fread(paste0(base_dir, "/overlap_hg19/cas9_overlap_matrix_gene.tsv"), header=TRUE, sep="\t")

barcode1_1 <- fread(paste0(base_dir, "/DMR_Results_hg19/hg19_barcode01.sorted.bed.gz_Liu_Jaenisch_FXS_supplemental1.bed.bed_stats.tsv"), header=TRUE, sep="\t")
barcode4_1 <- fread(paste0(base_dir, "/DMR_Results_hg19/hg19_barcode04.sorted.bed.gz_Liu_Jaenisch_FXS_supplemental1.bed.bed_stats.tsv"), header=TRUE, sep="\t")
barcode1_2 <- fread(paste0(base_dir, "/DMR_Results_hg19/EXP-NBD104_barcode01.bed.gz_Liu_Jaenisch_FXS_supplemental1.bed.bed_stats.tsv"), header=TRUE, sep="\t")
barcode4_2 <- fread(paste0(base_dir, "/DMR_Results_hg19/EXP-NBD104_barcode04.bed.gz_Liu_Jaenisch_FXS_supplemental1.bed.bed_stats.tsv"), header=TRUE, sep="\t")

#Let's make a new dataframe that includes chrom, start, end
combined <- barcode1_1[,c("chrom", "start", "end")]

#let's insert percent_m column from each barcode, with a column label
combined$barcode1_1_pm <- barcode1_1$percent_m
combined$barcode1_2_pm <- barcode1_2$percent_m
combined$barcode4_1_pm <- barcode4_1$percent_m
combined$barcode4_2_pm <- barcode4_2$percent_m
combined$barcode1_1_valid_m <- barcode1_1$count_valid_m
combined$barcode1_2_valid_m <- barcode1_2$count_valid_m
combined$barcode4_1_valid_m <- barcode4_1$count_valid_m
combined$barcode4_2_valid_m <- barcode4_2$count_valid_m

combined <- combined %>%
  mutate(
    mean_barcode1 = (barcode1_1_pm + barcode1_2_pm) / 2,
    mean_barcode4 = (barcode4_1_pm + barcode4_2_pm) / 2,
    delta = mean_barcode1 - mean_barcode4
  )


#Let's grab differentially methylated regions, by absolute difference of 0.3
dmr <- combined[abs(combined$delta) > 25,]

#Drop any sample with 0 valid_m in any barcodes
dmr <- dmr[dmr$barcode1_1_valid_m > 0 & dmr$barcode1_2_valid_m > 0 & dmr$barcode4_1_valid_m > 0 & dmr$barcode4_2_valid_m > 0,]

#let's save this to a file
write.csv(dmr, file = file.path(base_dir, "DMR_Results_hg19/CAS9_DMR.csv"))

#Let's also just save the dmr file as a bed file
write.table(dmr[,c("chrom", "start", "end")], file = file.path(base_dir, "DMR_Results_hg19/CAS9_DMR.bed"), sep="\t", quote=FALSE, row.names=FALSE)


















