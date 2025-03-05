#!/bin/bash

# Set paths to your BED graph files and genome sizes file
BEDGRAPH_DIR="/Users/carternorton/Desktop/barcoded_bams_simplex/bedgraphs"
GENOME_SIZES="/Users/carternorton/Desktop/hs1.chrom.sizes.txt"
OUTPUT_DIR="/Users/carternorton/Desktop/barcoded_bams_simplex/sorted_bedgraphs"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all BED graph files
for file in "$BEDGRAPH_DIR"/*.bed; do
  filename=$(basename "$file" .bed) # Extract the base name without extension

  # Sort the BED file and save it to the output directory
  sort -k1,1 -k2,2n "$file" > "${OUTPUT_DIR}/${filename}.sorted.bedgraph"

  # Compress the sorted BED graph file and create the index
  bgzip -c "${OUTPUT_DIR}/${filename}.sorted.bedgraph" > "${OUTPUT_DIR}/${filename}.bedgraph.gz"
  tabix -s 1 -b 2 -e 3 "${OUTPUT_DIR}/${filename}.bedgraph.gz"

  # Convert to TDF format using igvtools
  igvtools toTDF "${OUTPUT_DIR}/${filename}.sorted.bedgraph" "${OUTPUT_DIR}/${filename}.bedgraph.idf" "/Users/carternorton/Desktop/cluster_drive/ref/T2T-CHM13v2.fna"

  echo "Converted: ${filename}.bed"
done

echo "All files have been converted!"


