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
library(eulerr)

base_dir<- "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore"

#Let's read the gene table from overlaps
cas9_genes_all <- fread(paste0(base_dir, "/overlap_hg19/cas9_overlap_matrix_gene.tsv"), header=TRUE, sep="\t")
de_CGGFXS_vs_FXS <- fread(paste0(base_dir, "/RNA_HG01_hg19/FXS_CGG_ONLY_vs_dCAS9_CGG/hisat_deseq2_results.csv"), header=TRUE, sep=",")
de_FXS_vs_WT <- fread(paste0(base_dir, "/RNA_HG01_hg19/FXS_CGG_ONLY_vs_WT_dCAS9/hisat_deseq2_results.csv"), header=TRUE, sep=",")
de_NHG3FXS_vs_FXS <- fread(paste0(base_dir, "/RNA_HG01_hg19/FXS_CGG_ONLY_vs_NHG3_dCAS9/hisat_deseq2_results.csv"), header=TRUE, sep=",")
dmr_cas9_genes <- fread(paste0(base_dir, "/overlap_hg19/dmr_liu_jaenisch/overlap_matrix_gene.tsv"), header=TRUE, sep="\t")

#Rename the first column of both to gene_id
setnames(cas9_genes_all, "Gene", "gene")
setnames(dmr_cas9_genes, "Gene", "gene")
setnames(de_CGGFXS_vs_FXS, "V1", "gene")
setnames(de_NHG3FXS_vs_FXS, "V1", "gene")
setnames(de_FXS_vs_WT, "V1", "gene")

#Let's filter each to only include genes that are in each
intersect <- intersect(cas9_genes_all$gene, de_CGGFXS_vs_FXS$gene)
cas9_genes_all <- cas9_genes_all[cas9_genes_all$gene %in% intersect]
de_CGGFXS_vs_FXS <- de_CGGFXS_vs_FXS[de_CGGFXS_vs_FXS$gene %in% intersect]
print(length(intersect))

#Let's filter de_CGGFXS_vs_FXS to only include genes with a p-value < 0.001 and a log2FoldChange > 1.5
de_CGGFXS_vs_FXS_SIG<- de_CGGFXS_vs_FXS[de_CGGFXS_vs_FXS$padj < 0.05 & de_CGGFXS_vs_FXS$log2FoldChange >= 1,]
de_FXS_vs_WT_SIG <- de_FXS_vs_WT[de_FXS_vs_WT$padj < 0.05 & de_FXS_vs_WT$log2FoldChange >= 1,]
de_NHG3FXS_vs_FXS_SIG <- de_NHG3FXS_vs_FXS[de_NHG3FXS_vs_FXS$padj < 0.05 & de_NHG3FXS_vs_FXS$log2FoldChange >= 1,]

#Let's filter cas9_genes to only include genes with "ANY" > 0
cas9_genes <- cas9_genes_all[ANY > 0,]
dmr_cas9_genes <- dmr_cas9_genes[ANY > 0,]

#Let's get the length of de_CGGFXS_vs_FXS and length of cas9_genes
length(de_CGGFXS_vs_FXS_SIG$gene)
length(cas9_genes$gene)

#Let's get the intersection of the two gene lists
intersect_genes <- intersect(de_CGGFXS_vs_FXS_SIG$gene, cas9_genes$gene)
print(length(intersect_genes))
print(intersect_genes)

#Let's see where the Cas9 binding is on each of these
intersect_cas9 <- cas9_genes_all[cas9_genes_all$gene %in% intersect_genes,]
intersect_cas9_sums <- intersect_cas9[, 2:ncol(intersect_cas9)] %>% colSums()

#Let's make this plot
# Define the sets
set1 <- de_CGGFXS_vs_FXS_SIG$gene
set2 <- cas9_genes$gene

pdf(file = file.path(base_dir, "overlap_hg19/venn_plot.pdf")) 

venn.plot <- draw.pairwise.venn(
  area1 = length(set1),
  area2 = length(set2),
  cross.area = length(intersect(set1, set2)),
  category = c("DE Genes (padj < 0.05, log2fc > 1)", "Genes with Anti-Cas9 ChIP-Seq Peaks"),
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.pos = c(180, 0),  
  cat.dist = c(-0.05, -0.05)
)

grid.draw(venn.plot)
dev.off()

#Let's save the intersect genes to a file
write.table(intersect_genes, file=paste0(base_dir, "/overlap_hg19/cas9_DEFXSvFXSCAS9_intersect_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

#Let's see how many of these genes are TFs
# From https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt
tf_genes <- fread(paste0(base_dir, "/hg19_regions_of_interest/TFs.csv"), header=TRUE, sep=",", data.table=FALSE)
tf_genes <- tf_genes$"HGNC symbol"
intersect_tf_genes <- intersect(intersect_genes, tf_genes)

##Let's load in potential targets of these TFs, and see if any of them are in de non-intersect
# These are from https://chip-atlas.dbcls.jp, and include chip peaks within 5kb
runx2 <- fread(paste0(base_dir, "/hg19_regions_of_interest/RUNX2.tsv"), header=TRUE, sep="\t", data.table=FALSE)
ZNF703 <- fread(paste0(base_dir, "/hg19_regions_of_interest/ZNF703.tsv"), header=TRUE, sep="\t", data.table=FALSE)
RUNX2_genes <- runx2[runx2$`RUNX2|Average` > 250,]$"Target_genes"
ZNF703_genes <- ZNF703[ZNF703$`ZNF703|Average` > 250,]$"Target_genes"

print(length(RUNX2_genes))
print(length(ZNF703_genes))

#Let's get the intersection of non-intersect de genes and RUNX2_genes and ZNF703_genes
tf_targets <- unique(c(RUNX2_genes, ZNF703_genes))
de_only <- de_CGGFXS_vs_FXS_SIG[!de_CGGFXS_vs_FXS_SIG$gene %in% intersect_genes,]$gene
de_tf_targets <- intersect(tf_targets, de_only)
print(length(de_tf_targets))
#There are only 29 possible targets of TFs according to available chip-seq data
                     
### Now, let's read in PAR-CHIP and RIP-CHIP data
regions_folder <- paste0(base_dir, "/hg19_regions_of_interest/")
rip_clip <- fread(paste0(regions_folder, "Ascano_RIPCLIP-6.csv"), header=TRUE, sep=",", data.table=FALSE)
par_clip <- fread(paste0(regions_folder, "Ascano_PARCLIP-4.csv"), header=TRUE, sep=",", data.table=FALSE)

#Rename Columns
setnames(rip_clip, "RIP LFE (log2 fold enrichment)", "RIP_LF2")
setnames(rip_clip, "Supplementary Table 6 FLAG-HA FMRP iso1 RIP-chip Gene symbol", "gene")
setnames(par_clip, "FMR1 iso1_Reads", "iso1_reads")
setnames(par_clip, "FMR1 iso7_Reads", "iso7_reads")
setnames(par_clip, "FMR1 iso1_Total mRNA binding sites", "iso1_mRNA_binding_sites")
setnames(par_clip, "FMR1 iso7_Total mRNA binding sites", "iso7_mRNA_binding_sites")
setnames(par_clip, "Gene", "gene")

#Remove all genes with "NA" in the RIP_LF2 column, and over 2 LFC
rip_clip <- rip_clip[!is.na(rip_clip$RIP_LF2),]
rip_clip_genes <- rip_clip[rip_clip$RIP_LF2 >= 1,]$gene
#After finding the intersection of ripclip and par_clip (940), 
# only 19 genes were in the DE genes that were not in ANTI-CAS9 binding site.
# We will now stick to par_clip as a more lenient filter

#Let's get distribution of iso1_reads and iso1_mRNA_binding_sites
par_clip$total_binding <- par_clip$iso1_mRNA_binding_sites + par_clip$iso7_mRNA_binding_sites
print(summary(par_clip$total_binding))

#Let's select only genes with total_binding >= 2
par_clip_genes <- par_clip[par_clip$total_binding >= 2,]$gene

#Let's get the intersection of intersect_par_rip and de non-intersect_genes
de_only <- de_CGGFXS_vs_FXS_SIG[!de_CGGFXS_vs_FXS_SIG$gene %in% intersect_genes,]$gene
intersect_par_rip_gene<- intersect(par_clip_genes, de_only)
print(length(intersect_par_rip_gene))
#Only looks like 129 genes are differentialy expressed and have a binding site of FMRP-PAR-CLIP


#Let's see which of the intersect are also in the dmr genes
intersect_dmr_cas9_genes <- intersect(intersect_genes, dmr_cas9_genes$Gene)
print(length(intersect_dmr_cas9_genes))

#Let's get overlap between de_CGGFXS_vs_FXS_SIG and de_fxswt
intersect_fxswt_genes <- intersect(de_CGGFXS_vs_FXS_SIG$gene, de_FXS_vs_WT_SIG$gene)
print(length(intersect_fxswt_genes))

intersect_de_only_fxswt_genes <- intersect(de_only, de_FXS_vs_WT_SIG$gene)
print(length(intersect_de_only_fxswt_genes))

#Let's make a venn diagram of the overlap between de_CGGFXS_vs_FXS_SIG, de_FXS_vs_WT_SIG, and de_NHG3FXS_vs_FXS_SIG
pdf(file = file.path(base_dir, "overlap_hg19/venn_plot3.pdf"))

# Example sets
set1 <- de_CGGFXS_vs_FXS_SIG$gene
set2 <- de_FXS_vs_WT_SIG$gene
set3 <- de_NHG3FXS_vs_FXS_SIG$gene

# Compute overlaps
n12 <- length(intersect(set1, set2))
n13 <- length(intersect(set1, set3))
n23 <- length(intersect(set2, set3))
n123 <- length(Reduce(intersect, list(set1, set2, set3)))

venn.plot <- draw.triple.venn(
  area1 = length(set1),
  area2 = length(set2),
  area3 = length(set3),
  n12 = n12,
  n13 = n13,
  n23 = n23,
  n123 = n123,
  category = c("CGGFXS_vs_FXS_SIG", 
               "FXS_vs_WT_SIG",
               "NHG3FXS_vs_FXS_SIG"),
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cat.pos = c(180, 0, 270),  
  scaled = TRUE
)

grid.draw(venn.plot)
dev.off()


venn_data <- euler(c(
  "CGGFXS_vs_FXS" = length(set1),
  "WT_vs_FXS" = length(set2),
  "NHG3FXS_vs_FXS" = length(set3),
  "CGGFXS_vs_FXS&WT_vs_FXS" = n12,
  "CGGFXS_vs_FXS&NHG3FXS_vs_FXS" = n13,
  "WT_vs_FXS&NHG3FXS_vs_FXS" = n23,
  "CGGFXS_vs_FXS&WT_vs_FXS&NHG3FXS_vs_FXS" = n123
))

# Plot
pdf(file = file.path(base_dir, "overlap_hg19/venn_plot3_euler.pdf"))
plot(venn_data, fills = c("red", "blue", "green"), alpha = 0.5, quantities = TRUE)
dev.off()


### How many de genes are in NHG3
print(length(de_NHG3FXS_vs_FXS_SIG$gene))

#How many overlap with de_FXS_vs_WT
intersect_fxsnhg3_genes <- intersect(de_FXS_vs_WT_SIG$gene, de_NHG3FXS_vs_FXS_SIG$gene)
print(length(intersect_fxsnhg3_genes))








