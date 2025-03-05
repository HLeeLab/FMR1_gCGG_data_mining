if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
#BiocManager::install(c("tximport", "ggplot2", "DESeq2", "rtracklayer", "biomaRt")) # Add any other packages you need

#Lib path
lib_path <- "/home/cnorton5/R/x86_64-pc-linux-gnu-library/4.3"

#To install biomaRt you need to load system library 'libpng' on the cluster, then install R package png 
BiocManager::install(c("biomaRt"), lib = lib_path, dependencies = TRUE)

# Install CRAN packages if needed
#install.packages(c("tidyverse", "data.table"))

#Load Libraries
library(tximport)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(data.table)
library(rtracklayer)

base_dir = "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/hg19/rna_seq"
ref_dir = "/home/cnorton5/data_hlee308/cnorton5/ref"

#First, we need to generate a tx2gene object that maps transcript IDs to gene IDs
gtf <- import(paste(ref_dir, "hg19.67.ensGene.gtf", sep = "/"))
tx2gene <- as.data.frame(gtf)[, c("transcript_id", "gene_name")]
tx2gene <- unique(tx2gene)  # Remove duplicates
write.csv(tx2gene, file = paste(ref_dir, "hg19_other/tx2gene_hg19.67.csv", sep="/"),
                                row.names = FALSE)

# List of paths to quant.sf files
samples <- c("HG01mRNA-CGGonly_24d_mT_1/quant.sf", "HG01mRNA-CGGonly_24d_mT_2/quant.sf", 
             "HG01mRNA-dCas9E_CGG_1/quant.sf", "HG01mRNA-dCas9E_CGG_2/quant.sf")

#Append samples with base_dir/salmon_results
samples <- paste(base_dir,"salmon_results", samples, sep = "/")

# Import the quantification data
txi <- tximport(samples, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

coldata <- data.frame(
  row.names = samples,
  condition = c("control", "control", "treated", "treated")  # Adjust according to experiment
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
results <- results(dds)
results <- results[order(results$padj, na.last = NA), ]

# Save results to a CSV file
write.csv(as.data.frame(results), file = paste(base_dir,"deseq2_results.csv", sep="/"))

#Save txi
write.csv(as.data.frame(txi$counts), file = paste(base_dir,"txi_counts.csv", sep="/"))


## Plot PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pdf(paste(base_dir, "PCA_plot.pdf", sep = "/"))
ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")
dev.off()

##VOLCANO PLOT
volcano_data <- as.data.frame(results)
volcano_data$significance <- volcano_data$padj < 0.000001

if ("gene" %in% colnames(volcano_data)) {
  gene_col <- volcano_data$gene
} else {
  gene_col <- rownames(volcano_data)
}

# Add labels for points where -log10(padj) is above 20
volcano_data$label <- ifelse(-log10(volcano_data$padj) > 20 | gene_col == "FMR1", gene_col, "")

pdf(paste(base_dir, "volcano_plot.pdf", sep = "/"))
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point() +
  geom_text(aes(label = label), vjust = -0.5, hjust = 0.5, size = 3, check_overlap = TRUE) +  # Add labels
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") 
dev.off()

