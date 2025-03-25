# ---- Setup and Library Loading ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_pkgs <- c("VennDiagram", "Rsubread", "DESeq2", "ggplot2", "tidyverse",
                   "data.table", "dplyr", "eulerr")
suppressMessages(lapply(required_pkgs, library, character.only = TRUE))

base_dir <- "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore"

### HELPER FUNCTIONS ###
filter_by_sig <- function(dt, padj_cutoff = 0.05, lfc_cutoff = 1, direction = "up") {
  if (direction == "up") {
    dt[padj < padj_cutoff & log2FoldChange >= lfc_cutoff, ]
  } else {
    dt[padj < padj_cutoff & log2FoldChange <= -lfc_cutoff, ]
  }
}

# ---- Plotting Venn Diagrams ----
plot_venn_pair <- function(set1, set2, cat_labels, file_name) {
  pdf(file = file.path(base_dir, file_name))
  venn.plot <- draw.pairwise.venn(
    area1 = length(set1),
    area2 = length(set2),
    cross.area = length(intersect(set1, set2)),
    category = cat_labels,
    fill = c("red", "blue"),
    alpha = 0.5,
    cat.pos = c(180, 0),
    cat.dist = c(-0.05, -0.05)
  )
  grid.draw(venn.plot)
  dev.off()
}

get_euler <- function(set1, set2, set3, file_name) {
  
  # Euler diagram for triple DE comparisons
  # Compute unique counts (genes only in one category)
  n12 <- length(intersect(set1, set2))
  n13 <- length(intersect(set1, set3))
  n23 <- length(intersect(set2, set3))
  n123 <- length(Reduce(intersect, list(set1, set2, set3)))

  only_set1 <- length(set1) - (n12 + n13 - n123)
  only_set2 <- length(set2) - (n12 + n23 - n123)
  only_set3 <- length(set3) - (n13 + n23 - n123)

  # Compute Euler input with adjusted values
  venn_data <- euler(c(
    "CGGFXS_vs_FXS" = only_set1,
    "WT_vs_FXS" = only_set2,
    "NHG3FXS_vs_FXS" = only_set3,
    "CGGFXS_vs_FXS&WT_vs_FXS" = n12 - n123,  # Exclude triple overlap
    "CGGFXS_vs_FXS&NHG3FXS_vs_FXS" = n13 - n123,
    "WT_vs_FXS&NHG3FXS_vs_FXS" = n23 - n123,
    "CGGFXS_vs_FXS&WT_vs_FXS&NHG3FXS_vs_FXS" = n123
  ))

  # Save plot
  #pdf(file = file.path(base_dir, file_name))
  #plot(venn_data, fills = c("red", "blue", "green"), alpha = 0.5, quantities = TRUE)
  #dev.off()
  return(venn_data)
}


### END OF HELPER FUNCTIONS ###

# ---- Data Loading ----
cas9_genes_all <- fread(file.path(base_dir, "overlap_hg19/cas9_overlap_matrix_gene.tsv"), header = TRUE, sep = "\t")
dmr_cas9_genes <- fread(file.path(base_dir, "overlap_hg19/dmr_liu_jaenisch/overlap_matrix_gene.tsv"), header = TRUE, sep = "\t")
de_CGGFXS_vs_FXS <- fread(file.path(base_dir, "RNA_HG01_hg19/FXS_CGG_ONLY_vs_dCAS9_CGG/hisat_deseq2_results.csv"), header = TRUE, sep = ",")
de_FXS_vs_WT <- fread(file.path(base_dir, "RNA_HG01_hg19/FXS_CGG_ONLY_vs_WT_dCAS9/hisat_deseq2_results.csv"), header = TRUE, sep = ",")
de_NHG3FXS_vs_FXS <- fread(file.path(base_dir, "RNA_HG01_hg19/FXS_CGG_ONLY_vs_NHG3_dCAS9/hisat_deseq2_results.csv"), header = TRUE, sep = ",")

# Rename gene columns for consistency
setnames(cas9_genes_all, "Gene", "gene")
setnames(dmr_cas9_genes, "Gene", "gene")
setnames(de_CGGFXS_vs_FXS, "V1", "gene")
setnames(de_FXS_vs_WT, "V1", "gene")
setnames(de_NHG3FXS_vs_FXS, "V1", "gene")

# ---- Filtering and Intersections ----
# Filter to common genes between cas9 binding and DE data
common_genes <- intersect(cas9_genes_all$gene, de_CGGFXS_vs_FXS$gene)
cas9_genes_all <- cas9_genes_all[gene %in% common_genes]
de_CGGFXS_vs_FXS <- de_CGGFXS_vs_FXS[gene %in% common_genes]

# Filter DE genes (up-regulated)
de_CGGFXS_vs_FXS_SIG <- filter_by_sig(de_CGGFXS_vs_FXS)
de_FXS_vs_WT_SIG     <- filter_by_sig(de_FXS_vs_WT)
de_NHG3FXS_vs_FXS_SIG <- filter_by_sig(de_NHG3FXS_vs_FXS)

# Further filter cas9 and DMR genes for ANY > 0
cas9_genes <- cas9_genes_all[ANY > 0, ]
dmr_cas9_genes <- dmr_cas9_genes[ANY > 0, ]

# Intersection between DE genes and cas9 binding genes
intersect_genes <- intersect(de_CGGFXS_vs_FXS_SIG$gene, cas9_genes$gene)
cat("Intersected gene count:", length(intersect_genes), "\n")
print(intersect_genes)

# Sum Cas9 binding signals for intersected genes (excluding gene column)
intersect_cas9 <- cas9_genes_all[gene %in% intersect_genes, ]
intersect_cas9_sums <- colSums(intersect_cas9[, -1])
print(intersect_cas9_sums)

# Save intersected genes to a file
write.table(intersect_genes, file = file.path(base_dir, "overlap_hg19/cas9_DEFXSvFXSCAS9_intersect_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- TF Analysis ----
tf_genes <- fread(file.path(base_dir, "hg19_regions_of_interest/TFs.csv"), header = TRUE, sep = ",", data.table = FALSE)[, "HGNC symbol"]
intersect_tf_genes <- intersect(intersect_genes, tf_genes)

# Read potential TF target data
runx2 <- fread(file.path(base_dir, "hg19_regions_of_interest/RUNX2.tsv"), header = TRUE, sep = "\t", data.table = FALSE)
ZNF703 <- fread(file.path(base_dir, "hg19_regions_of_interest/ZNF703.tsv"), header = TRUE, sep = "\t", data.table = FALSE)
RUNX2_genes <- runx2[runx2$`RUNX2|Average` > 250, "Target_genes"]
ZNF703_genes <- ZNF703[ZNF703$`ZNF703|Average` > 250, "Target_genes"]

cat("RUNX2 targets:", length(RUNX2_genes), "\n")
cat("ZNF703 targets:", length(ZNF703_genes), "\n")

tf_targets <- unique(c(RUNX2_genes, ZNF703_genes))
de_only <- setdiff(de_CGGFXS_vs_FXS_SIG$gene, intersect_genes)
de_tf_targets <- intersect(tf_targets, de_only)
cat("DE TF target count:", length(de_tf_targets), "\n")

# ---- PAR-CLIP and RIP-CLIP Data ----
regions_folder <- file.path(base_dir, "hg19_regions_of_interest/")
rip_clip <- fread(file.path(regions_folder, "Ascano_RIPCLIP-6.csv"), header = TRUE, sep = ",", data.table = FALSE)
par_clip <- fread(file.path(regions_folder, "Ascano_PARCLIP-4.csv"), header = TRUE, sep = ",", data.table = FALSE)

# Rename columns for clarity
setnames(rip_clip, "Supplementary Table 6 FLAG-HA FMRP iso1 RIP-chip Gene symbol", "gene")
setnames(rip_clip, "RIP LFE (log2 fold enrichment)", "RIP_LF2")
setnames(par_clip, "Gene", "gene")
setnames(par_clip,
         c("FMR1 iso1_Reads", "FMR1 iso7_Reads",
           "FMR1 iso1_Total mRNA binding sites", "FMR1 iso7_Total mRNA binding sites"),
         c("iso1_reads", "iso7_reads", "iso1_mRNA_binding_sites", "iso7_mRNA_binding_sites"))

# Filter RIP-CLIP data and calculate PAR-CLIP binding totals
rip_clip <- rip_clip[!is.na(RIP_LF2) & RIP_LF2 >= 1, ]
par_clip$total_binding <- par_clip$iso1_mRNA_binding_sites + par_clip$iso7_mRNA_binding_sites
cat("Total binding summary:\n")
print(summary(par_clip$total_binding))
par_clip_genes <- par_clip[par_clip$total_binding >= 6, "gene"]

# Intersection of PAR/RIP binding genes with DE genes not in Cas9 peaks
de_only <- setdiff(de_CGGFXS_vs_FXS_SIG$gene, intersect_genes)
intersect_par_rip_gene <- intersect(par_clip_genes, de_only)
cat("Intersection of PAR/RIP genes with DE only:", length(intersect_par_rip_gene), "\n")

# ---- DMR and Additional DE Overlaps ----
intersect_dmr_cas9_genes <- intersect(intersect_genes, dmr_cas9_genes$gene)
cat("DMR & Cas9 intersect count:", length(intersect_dmr_cas9_genes), "\n")

intersect_fxswt_genes <- intersect(de_CGGFXS_vs_FXS_SIG$gene, de_FXS_vs_WT_SIG$gene)
cat("Overlap between CGGFXS_vs_FXS_SIG & FXS_vs_WT_SIG:", length(intersect_fxswt_genes), "\n")

intersect_de_only_fxswt_genes <- intersect(de_only, de_FXS_vs_WT_SIG$gene)
cat("Overlap between DE-only and FXS_vs_WT_SIG:", length(intersect_de_only_fxswt_genes), "\n")

### ---- Up-Regulated Genes ----
set1 <- de_CGGFXS_vs_FXS_SIG$gene
set2 <- de_FXS_vs_WT_SIG$gene
set3 <- de_NHG3FXS_vs_FXS_SIG$gene

#Euler
euler_data <- get_euler(set1, set2, set3, "overlap_hg19/venn_plot3_euler_up.pdf")
pdf(file = file.path(base_dir, "overlap_hg19/venn_plot3_euler_up.pdf"))
plot(euler_data, fills = c("red", "blue", "green"), alpha = 0.5, quantities = TRUE)
dev.off()

# Venn diagram for DE genes vs. Cas9 binding genes (up-regulated)
set_DE <- de_CGGFXS_vs_FXS_SIG$gene
set_Cas9 <- cas9_genes$gene
plot_venn_pair(set_DE, set_Cas9,
               c("DE Genes (padj < 0.05, log2fc > 1)", "Genes with Anti-Cas9 ChIP-Seq Peaks"),
               "overlap_hg19/venn_plot_up_CAS9_bind.pdf")


### ---- Down-Regulated Genes ----
de_CGGFXS_vs_FXS_SIG_Down <- filter_by_sig(de_CGGFXS_vs_FXS, direction = "down")
de_FXS_vs_WT_SIG_Down     <- filter_by_sig(de_FXS_vs_WT, direction = "down")
de_NHG3FXS_vs_FXS_SIG_Down <- filter_by_sig(de_NHG3FXS_vs_FXS, direction = "down")

set1_down <- de_CGGFXS_vs_FXS_SIG_Down$gene
set2_down <- de_FXS_vs_WT_SIG_Down$gene
set3_down <- de_NHG3FXS_vs_FXS_SIG_Down$gene

#Euler
down_euler_data <- get_euler(set1_down, set2_down, set3_down, "overlap_hg19/venn_plot3_euler_down.pdf")
pdf(file = file.path(base_dir, "overlap_hg19/venn_plot3_euler_down.pdf"))
plot(down_euler_data, fills = c("red", "blue", "green"), alpha = 0.5, quantities = TRUE)
dev.off()

#Venn Pair
plot_venn_pair(set1_down, set_Cas9,
               c("DE Genes Down (padj < 0.05, log2fc <= -1)", "Genes with Anti-Cas9 ChIP-Seq Peaks"),
               "overlap_hg19/venn_plot_down_CAS9_bind.pdf")

de_only_down <- setdiff(de_CGGFXS_vs_FXS_SIG_Down$gene, intersect_genes)
intersect_par_rip_gene_down <- intersect(par_clip_genes, de_only_down)
cat("Intersection of PAR/RIP genes with DE-only down:", length(intersect_par_rip_gene_down), "\n")







