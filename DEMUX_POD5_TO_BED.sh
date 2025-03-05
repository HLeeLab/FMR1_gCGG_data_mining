#!/bin/bash

#SBATCH --partition=l40s         # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1                # Number of nodes (1 node)
#SBATCH --ntasks-per-node=48     # Number of tasks per node (48 cores)
#SBATCH --mem=96G                # Memory allocation (96 GB)
#SBATCH --gres=gpu:4             # Request 4 GPUs
#SBATCH -A hlee308_gpu           # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/DEMUX_POD5_To_BAM.log  # Standard output log file (with job ID)
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/DEMUX_POD5_To_BAM_error.log  # Error log file

# Set the working directory 
cd /home/cnorton5 || exit 1

# Base directories
base_dir=/home/cnorton5/scr4_hlee308
project_dir="$base_dir/cnorton5/old_nanopore/ONT_HG01_hg19"
pod5_dir="$base_dir/cnorton5/old_nanopore/ONT_HG01_hg19/pod5"
ref_file=/home/cnorton5/data_hlee308/cnorton5/ref/hg19.fa

module load samtools

# Check for executables
command -v pod5 >/dev/null 2>&1 || { echo "pod5 not found"; exit 1; }
command -v dorado >/dev/null 2>&1 || { echo "dorado not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v modkit >/dev/null 2>&1 || { echo "modkit not found"; exit 1; }

# Check for reference file
[ -f "$ref_file" ] || { echo "Reference file not found: $ref_file"; exit 1; }

# Make necessary directories
mkdir -p "$project_dir"
mkdir -p "$pod5_dir"

#Given a barcode directory, find all subdirectories and process each one
input_dir="$base_dir/ONT-HG01_fast5_DN"

# Find all subdirectories that contain "barcode" in their name
barcode_dirs=($(find "$input_dir" -type d -name "*barcode*"))

#To test, print out the barcode directories
echo "Barcode directories:"
for dir in "${barcode_dirs[@]}"; do
  echo "$dir"
done

# Loop over each barcode directory
for input_dir in "${barcode_dirs[@]}"; do
  echo "=============================="
  echo "Processing directory: $input_dir"

  barcode_name=$(basename "$input_dir")
  
  # Check if the input directory exists
  if [ ! -d "$input_dir" ]; then
    echo "Input directory does not exist: $input_dir. Skipping..."
    continue
  fi

  # Enable nullglob so non-matching globs are removed
  shopt -s nullglob
  pod5_files=("$input_dir"/*.pod5)
  fast5_files=("$input_dir"/*.fast5)
  merged_pod5="$pod5_dir/${barcode_name}_merged.pod5"
  
  # Merge or convert files if necessary
  if [ -f "$merged_pod5" ]; then
    echo "Found merged.pod5 in $input_dir ... continuing with alignment"
  elif [ ${#pod5_files[@]} -gt 0 ]; then
    echo "Found pod5 files in $input_dir ... merging into merged.pod5"
    pod5 merge "${pod5_files[@]}" --output "$merged_pod5" --threads 48
  elif [ ${#fast5_files[@]} -gt 0 ]; then
    echo "Found fast5 files in $input_dir ... converting to pod5 and merging into merged.pod5"
    pod5 convert fast5 "${fast5_files[@]}" --output "$merged_pod5" --threads 48
  else 
    echo "No fast5 or pod5 files found in $input_dir. Skipping..."
    continue
  fi
  
  # Determine the output BAM filename (includes barcode directory name)
  align_genome=$(basename "$ref_file" | awk -F. '{print $1}')
  output_name="${align_genome}_${barcode_name}.bam"
  output_path="$project_dir/$output_name"
  
  # Run Dorado basecaller with GPUs
  echo "Running dorado basecaller for $input_dir..."
  dorado basecaller \
    /home/cnorton5/data_hlee308/cnorton5/usr/bin/dorado-0.9.0-linux-x64/bin/models/dna_r9.4.1_e8_hac@v3.3 \
    "$merged_pod5" \
    --modified-bases 5mCG_5hmCG \
    --kit-name EXP-NBD104 \
    --reference "$ref_file" > "$output_path"
  
  # Since data are already demuxed, we skip the demux step
  
  # Process the generated BAM file
  if [ -f "$output_path" ]; then
    echo "Processing BAM file: $output_path"
    sorted_bam="${output_path%.bam}.sorted.bam"
    
    # Sort the BAM file
    samtools sort --threads 24 -o "$sorted_bam" "$output_path"
    rm "$output_path"
    
    # Index the sorted BAM file
    samtools index "$sorted_bam" -@ 24
    
    # Run modkit pileup
    modkit pileup "$sorted_bam" "${sorted_bam%.bam}.bed" --threads 24

  else
    echo "BAM file not found: $output_path. Skipping processing for this directory."
  fi
  
  echo "Finished processing directory: $input_dir"
done

echo "=============================="
echo "All processing complete!"

