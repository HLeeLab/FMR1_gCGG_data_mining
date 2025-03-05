#!/bin/bash

# Set paths to your BED graph files and genome sizes file
BEDGRAPH_DIR="/Users/carternorton/Desktop/barcoded_bams_simplex/bedgraphs"
GENOME_SIZES="/Users/carternorton/Desktop/hs1.chrom.sizes.txt"
OUTPUT_DIR="/Users/carternorton/Desktop/barcoded_bams_simplex/bigwigs"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all BED graph files
for file in "$BEDGRAPH_DIR"/*.bed; do
  # Sort the BED graph file
  sorted_file="${file%.bedgraph}_sorted.bedgraph"
  #Remove the 5th column (num of cpgs) and sort
  awk '{ $5=""; print $0 }' "$file" | awk '{$1=$1; print}' | sort -k1,1 -k2,2n > "$sorted_file"

  # Convert to BigWig
  output_bw="$OUTPUT_DIR/$(basename "${file%.bedgraph}.bw")"
  /Users/carternorton/Desktop/bedGraphToBigWig "$sorted_file" "$GENOME_SIZES" "$output_bw"

  # Optionally, delete the sorted intermediate file
  rm "$sorted_file"

  echo "Converted: $file -> $output_bw"
done

echo "All files have been converted!"

