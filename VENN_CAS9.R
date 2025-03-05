# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor and CRAN packages
#BiocManager::install(c("Rsubread")) #, "DESeq2", "ggplot2", "tidyverse", "data.table"), ask = FALSE)

# Load libraries
library(VennDiagram)
library(Rsubread)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(data.table)
library(dplyr)

base_dir<- "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore"

#Let's read the gene table from overlaps
cas9_genes_all <- fread(paste0(base_dir, "/overlap_hg19/overlap_matrix_gene.tsv"), header=TRUE, sep="\t")
de_genes_all <- fread(paste0(base_dir, "/RNA_HG01_hg19/hisat_deseq2_results.csv"), header=TRUE, sep=",")

#Rename the first column of both to gene_id
setnames(cas9_genes_all, "Gene", "gene")
setnames(de_genes_all, "V1", "gene")

#Let's filter de_genes to only include genes with a p-value < 0.001 and a log2FoldChange > 1.5
de_genes <- de_genes_all[de_genes_all$padj < 0.01 & de_genes_all$log2FoldChange > 1,]

#Let's filter cas9_genes to only include genes with "ANY" > 0
cas9_genes <- cas9_genes_all[ANY > 0,]

#Let's get the length of de_genes and length of cas9_genes
length(de_genes$gene)
length(cas9_genes$gene)

#Let's get the intersection of the two gene lists
intersect_genes <- intersect(de_genes$gene, cas9_genes$gene)

print(length(intersect_genes))

print(intersect_genes)

#How many of intersect genes have a 1 value in the "UP" column of cas9_genes?
intersect_cas9 <- cas9_genes_all[cas9_genes_all$gene %in% intersect_genes,]

#Let's make a bargraph of the sum of each column in intersect_cas9
intersect_cas9_sums <- intersect_cas9[, 2:ncol(intersect_cas9)] %>% colSums()


#Let's make this plot
# Define the sets
set1 <- de_genes$gene
set2 <- cas9_genes$gene

# Create the Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(set1),
  area2 = length(set2),
  cross.area = length(intersect(set1, set2)),
  category = c("DE Genes ( padj < 0.01, lfc > 1)", "Genes with Anti-Cas9 ChIP-Seq Peaks"),
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.pos = c(180, 0),  # Positions labels toward the center
  cat.dist = c(-0.05, -0.05)  # Moves labels inward
)

# Display the Venn diagram]
pdf(file = file.path(base_dir, "overlap_hg19/venn_plot.pdf")) 
grid.draw(venn.plot)
dev.off()




