# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor and CRAN packages
BiocManager::install(c("Rsubread")) #, "DESeq2", "ggplot2", "tidyverse", "data.table"), ask = FALSE)

# Load libraries
library(Rsubread)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(data.table)

# Define directories (adjust these paths to your setup)
base_dir <- "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/hg19/rna_seq"
bam_dir  <- file.path(base_dir, "hisat2_results")  # Directory containing BAM files
ref_dir  <- "/home/cnorton5/data_hlee308/cnorton5/ref"

# Define path to the GTF annotation file
annotation_file <- file.path(ref_dir, "hg19.67.ensGene.gtf")

# List BAM files (assumes files end with .bam)
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE, recursive=TRUE)
# Name the vector using the file base names (without directory path)
names(bam_files) <- basename(bam_files)

# Run featureCounts using Rsubread
# We set:
#   - isGTFAnnotationFile=TRUE to indicate a GTF file is provided
#   - GTF.featureType="exon" to count reads overlapping exons
#   - GTF.attrType="gene_name" to use the gene name attribute as the gene identifier
#   - useMetaFeatures=TRUE to collapse counts from all exons of the same gene
fc <- featureCounts(files = bam_files,
                    annot.ext = annotation_file,
                    isGTFAnnotationFile = TRUE,
                    GTF.featureType = "exon",
                    GTF.attrType = "gene_name",  # Change to "gene_id" if desired
                    useMetaFeatures = TRUE,
                    nthreads = 12, 
                    isPairedEnd=TRUE)  # Adjust number of threads as needed

# Extract the counts matrix
counts <- fc$counts

#Save the feature counts table
write.csv(as.data.frame(counts),
          file = file.path(base_dir, "/hisat2_counts.csv"))

#Read csv
counts <- fread(file = file.path(base_dir, "/hisat2_counts.csv"), header = TRUE, data.table = FALSE)

#Set row names to the first column
row.names(counts) <- counts$V1
counts <- counts[, -1]

# Create sample metadata
# Ensure that the order of samples in coldata matches the columns in 'counts'
# Here we assume that the BAM file names correspond to the sample names
sample_names <- names(bam_files)
coldata <- data.frame(
  row.names = sample_names,
  condition = c("control", "control", "treated", "treated")  # Adjust based on your experiment
)

# Create DESeq2 dataset from the counts matrix and metadata
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Run the differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj, na.last = NA), ]

# Save DESeq2 results and counts to CSV files
write.csv(as.data.frame(res),
          file = file.path(base_dir, "hisat_deseq2_results.csv"))

# --------------------
# PCA Plot
# --------------------
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pdf(file = file.path(base_dir, "hisat2_PCA_plot.pdf"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")
dev.off()

# --------------------
# Volcano Plot
# --------------------
volcano_data <- as.data.frame(res)
volcano_data$significance <- volcano_data$padj < 0.000001

# Use rownames as gene names (these were set from the 'gene_name' in the annotation)
volcano_data$gene <- rownames(volcano_data)
# Add labels for points with high significance or a specific gene of interest ("FMR1")
volcano_data$label <- ifelse(-log10(volcano_data$padj) > 20 | volcano_data$gene == "FMR1",
                             volcano_data$gene, "")

pdf(file = file.path(base_dir, "hisat2_volcano_plot.pdf"))
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point() +
  geom_text(aes(label = label), vjust = -0.5, hjust = 0.5, size = 3, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")
dev.off()

