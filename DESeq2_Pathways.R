# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor and CRAN packages
BiocManager::install(c("clusterProfiler", "msigdbr", "ggplot2", "dplyr"))
BiocManager::install(c("enrichplot"))
# Load required libraries
library(clusterProfiler)  # For GSEA
library(msigdbr)          # For downloading gene sets from MSigDB
library(ggplot2)          # For plotting
library(dplyr)            # For data manipulation
library(enrichplot)       # For plotting GSEA results

# Function to read DESeq2 results, perform GSEA, and generate a plot
performGSEA <- function(file, 
                        species = "Homo sapiens", 
                        gene_column = "gene", 
                        logfc_column = "log2FoldChange", 
                        pval_column = "pvalue",
                        gs_category = "C5",        # e.g., "H" for hallmark gene sets
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        pvalueCutoff = 0.05,
                        output_dir = "") {
  
  #Set output dir
  output_pdf <- paste0(output_dir, "GSEA_plot.pdf")
  output_GSEA <- paste0(output_dir, "GSEA_result.csv")
  
  # Read DESeq2 results (adjust the read function if your file is TSV or other)
  res <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  print(head(res))
  
  #rename the first column "gene"
  colnames(res)[colnames(res) == "X"] <- "gene"
  
  
  # Remove rows with missing log fold changes
  res <- res[!is.na(res[[logfc_column]]), ]
  
  # Create a named vector of log2 fold changes.
  # The names should be gene symbols (or IDs) as expected by the gene sets.
  geneList <- res[[logfc_column]]
  names(geneList) <- res[[gene_column]]
  
  # Sort in decreasing order (required by GSEA)
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Download gene sets for the specified species and category using msigdbr.
  msigdbr_df <- msigdbr(species = species)
  
  # Filter for the desired gene set category (e.g., hallmark "H")
  gene_sets_df <- msigdbr_df %>% filter(gs_cat == gs_category)
  
  # Convert to a list where each element is a gene set (a vector of gene symbols)
  gene_sets <- gene_sets_df %>% 
    split(x = .$gene_symbol, f = .$gs_name)
  
  # Prepare TERM2GENE data.frame for clusterProfiler::GSEA.
  # Here, we need two columns: one with the term (gene set name) and one with gene symbols.
  term2gene <- stack(gene_sets)[, c("ind", "values")]
  colnames(term2gene) <- c("term", "gene")
  
  # Perform GSEA
  gseaRes <- GSEA(geneList,
                  TERM2GENE = term2gene,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  pvalueCutoff = pvalueCutoff)
  
  # Convert the GSEA results to a data frame for plotting
  gsea_df <- as.data.frame(gseaRes)
  
  # For demonstration, we use the normalized enrichment score (NES) as a proxy for fold enrichment.
  # (Note: In ORA analyses, fold enrichment is often computed as observed/expected gene counts.)
  gsea_df$FoldEnrichment <- gsea_df$NES
  
  #save the GSEA result
  write.csv(gsea_df, file = output_GSEA, row.names = FALSE)
  
  # Return the GSEA results and the plot object
  return(gseaRes)
}

# Define directories (adjust these paths to your setup)
base_dir <- "/Users/carternorton/Desktop/cluster_sdrive/old_nanopore/RNA_HG01_hg19/"

# Example usage
# Replace "DESeq2_results.csv" with the path to your DESeq2 results file.

fxsCGG_vs_dCas9_path <- paste0(base_dir, "FXS_CGG_ONLY_vs_dCAS9_CGG/")
deseq_path <- paste0(fxsCGG_vs_dCas9_path, "hisat_deseq2_results.csv")
result <- performGSEA(deseq_path, output_dir = fxsCGG_vs_dCas9_path)

pdf(file=paste0(fxsCGG_vs_dCas9_path, "GSEA_plot.pdf"), width = 10, height = 6)
ridgeplot(result, showCategory=10, fill="p.adjust", core_enrichment = TRUE)
dev.off()

FXS_CGG_ONLY_vs_NHG3_dCAS9_path <- paste0(base_dir, "FXS_CGG_ONLY_vs_NHG3_dCAS9/")
deseq_path <- paste0(FXS_CGG_ONLY_vs_NHG3_dCAS9_path, "hisat_deseq2_results.csv")
result <- performGSEA(deseq_path, output_dir = FXS_CGG_ONLY_vs_NHG3_dCAS9_path)

pdf(file=paste0(FXS_CGG_ONLY_vs_NHG3_dCAS9_path, "GSEA_plot.pdf"), width = 10, height = 6)
ridgeplot(result, showCategory=10, fill="p.adjust", core_enrichment = TRUE)
dev.off()

FXS_CGG_ONLY_vs_WT_dCAS9_path <- paste0(base_dir, "FXS_CGG_ONLY_vs_WT_dCAS9/")
deseq_path <- paste0(FXS_CGG_ONLY_vs_WT_dCAS9_path, "hisat_deseq2_results.csv")
result <- performGSEA(deseq_path, output_dir = FXS_CGG_ONLY_vs_WT_dCAS9_path)

pdf(file=paste0(FXS_CGG_ONLY_vs_WT_dCAS9_path, "GSEA_plot.pdf"), width = 10, height = 6)
ridgeplot(result, showCategory=10, fill="p.adjust", core_enrichment = TRUE)
dev.off()




