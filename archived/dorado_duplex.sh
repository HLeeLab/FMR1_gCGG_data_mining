#!/bin/bash

#SBATCH --partition=l40s # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=48
#SBATCH --mem=96G # Memory allocation (96 GB)
#SBATCH --gres=gpu:4 # 
#SBATCH -A hlee308_gpu # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/dorado_output.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/dorado_error.log # Error log file

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1

# Set base directories
base_dir=/home/cnorton5/scr4_hlee308
nanopore_dir="$base_dir/cnorton5/old_nanopore"
pod5_dir="$nanopore_dir/full_pod5"
barcoded_bams_dir="$nanopore_dir/barcoded_bams_duplex"
methylbed_dir="$barcoded_bams_dir/methylbed"

module load samtools

#Check for executables
command -v pod5 >/dev/null 2>&1 || { echo "pod5 not found"; exit 1; }
command -v dorado >/dev/null 2>&1 || { echo "dorado not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v modkit >/dev/null 2>&1 || { echo "modkit not found"; exit 1; }
[ -f /home/cnorton5/data_hlee308/cnorton5/ref/T2T-CHM13v2.fna ] || { echo "Reference file not found: $base_dir/ref/T2T-CHM13v2.fna"; exit 1; }

# Make necessary directories
mkdir -p "$nanopore_dir"
mkdir -p "$pod5_dir"
mkdir -p "$barcoded_bams_dir"
mkdir -p "$methylbed_dir"

# # First, move and convert all fast5 files
# pod5 convert fast5 /home/cnorton5/data_hlee308/ONT-HG02_fast5_pass/*.fast5 --output "$pod5_dir" 
# #--one-to-one home/cnorton5/data_hlee308/ONT-HG02_fast5_pass/*.fast5

# # Create channels for POD5 data
# pod5 view "$pod5_dir" --include "read_id, channel" --output "$pod5_dir/summary.tsv"
# pod5 subset "$pod5_dir" --summary "$pod5_dir/summary.tsv" --columns channel --output "$pod5_dir/split_by_channel"

# Command to run Dorado basecaller with GPUs
dorado duplex \
  /home/cnorton5/data_hlee308/cnorton5/usr/bin/dorado-0.9.0-linux-x64/bin/models/dna_r9.4.1_e8_hac@v3.3 \
  "$pod5_dir/split_by_channel" \
  --modified-bases 5mCG_5hmCG \
  --kit-name EXP-NBD104 \
  --reference /home/cnorton5/data_hlee308/cnorton5/ref/T2T-CHM13v2.fna > "$nanopore_dir/full_output_duplex.bam"

# Separate by barcode
dorado demux --output-dir "$barcoded_bams_dir" --no-classify "$nanopore_dir/full_output_duplex.bam"

#rm "$nanopore_dir/full_output.bam"

# Process each BAM file
for bam_file in "$barcoded_bams_dir"/*.bam; do 
  # Define sorted BAM filename
  sorted_bam_file="${bam_file%.bam}.sorted.bam"

  # Sort the BAM file
  samtools sort -o "$sorted_bam_file" "$bam_file"

  rm "$bam_file"
  
  # Index the sorted BAM file
  samtools index "$sorted_bam_file"
  
  # Run modkit on the sorted BAM file
  modkit pileup "$sorted_bam_file" "${bam_file%.bam}.bed"
done

#Let's remove the prefix from all files in barcoded_bams_dir, removing everything before the first underscore
for file in "$barcoded_bams_dir"/*; do
  mv "$file" "$(echo "$file" | sed 's/^[^_]*_//')"
done

echo "Processing complete!"
