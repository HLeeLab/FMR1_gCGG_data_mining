# Load required packages
#BiocManager::install("Rsamtools")
#install.packages("R.utils")
#install.packages("data.table")
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(R.utils)
library(data.table)

# Define paths to chain files
chain_hg19 <- "/home/cnorton5/data_hlee308/cnorton5/ref/liftover/hs1ToHg19.over.chain"
chain_hg38 <- "/home/cnorton5/data_hlee308/cnorton5/ref/liftover/hs1ToHg38.over.chain"

# Load chain files
hg19_chain <- import.chain(chain_hg19)
hg38_chain <- import.chain(chain_hg38)

# Define input and output directories
base_dir <- "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/barcoded_bams_simplex"
input_dir <- "/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/barcoded_bams_simplex/sorted_bedgraphs"  # Update with your input folder

output_hg19_dir <- paste(base_dir, "output_hg19_dir", sep="/")
output_hg38_dir <- paste(base_dir, "output_hg38_dir", sep="/")

# Create output directories if they don't exist
dir.create(output_hg19_dir)
dir.create(output_hg38_dir)

# List of all .gz bedGraph files in the input directory
bedgraph_files <- list.files(input_dir, pattern = "\\.bedgraph$", full.names = TRUE)

# Load the chromosome mapping file
chromosome_mapping <- fread("/home/cnorton5/data_hlee308/cnorton5/ref/T2T_other/T2T_chrom_mapping.tsv",
                            sep = "\t")
# Ensure only the relevant columns are kept
chromosome_mapping <- chromosome_mapping[, c("RefSeq seq accession", "UCSC style name")]

# Create a named vector for mapping
mapping_vec <- setNames(chromosome_mapping$`UCSC style name`,
                        chromosome_mapping$`RefSeq seq accession`)

for (file in bedgraph_files) {
  
  # Read in the bedGraph file
  bed_data <- fread(file, sep= "\t")
  bed_data <- GRanges(seqnames = bed_data$V1, ranges = IRanges(start = bed_data$V2, end = bed_data$V3), score = bed_data$V4)
  
  # Convert seqnames
  old_seqnames <- as.character(seqnames(bed_data))
  new_seqnames <- mapping_vec[old_seqnames]
  new_seqnames[is.na(new_seqnames)] <- old_seqnames[is.na(new_seqnames)]
  
  # Assign new seqnames
  seqlevels(bed_data) <- union(seqlevels(bed_data), unique(na.omit(new_seqnames)))
  seqnames(bed_data) <- factor(new_seqnames, levels = seqlevels(bed_data))
  
  cat("Memory before liftover:", system("grep 'MemAvailable' /proc/meminfo", intern = TRUE), "\n", file = stderr())
  
  # Liftover to hg19
  lifted_hg19 <- liftOver(bed_data, hg19_chain)
  lifted_hg19 <- unlist(lifted_hg19)  # Convert list to GRanges object
  hg19_output_file <- file.path(output_hg19_dir, paste0(basename(file), "_hg19.bedGraph"))
  fwrite(as.data.table(lifted_hg19), file = hg19_output_file, sep = "\t", quote = FALSE, col.names = FALSE)
  
  # Liftover to hg38
  lifted_hg38 <- liftOver(bed_data, hg38_chain)
  lifted_hg38 <- unlist(lifted_hg38)  # Convert list to GRanges object
  hg38_output_file <- file.path(output_hg38_dir, paste0(basename(file), "_hg38.bedGraph"))
  fwrite(as.data.table(lifted_hg38), file = hg38_output_file, sep = "\t", quote = FALSE, col.names = FALSE)
  
  cat("Liftover completed for:", basename(file), "\n")
  # Force memory cleanup
  bed_data <- NULL
  lifted_hg19 <- NULL
  lifted_hg38 <- NULL
  new_seqnames <- NULL
  old_seqnames <- NULL
  gc()
  cat("Memory after liftover:", system("grep 'MemAvailable' /proc/meminfo", intern = TRUE), "\n", file = stderr())
}




