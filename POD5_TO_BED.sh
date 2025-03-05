#!/bin/bash

#SBATCH --partition=l40s # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=48 # Number of tasks per node (24 cores)
#SBATCH --gres=gpu:4 # 
#SBATCH -A hlee308_gpu # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/POD5_To_BAM.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/POD5_To_BAM_error.log # Error log file

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1

# Set base directories
base_dir="/home/cnorton5/scr4_hlee308/"
project_dir="$base_dir/cnorton5/old_nanopore/ONT_HG02_hg19_full"
input_dir="$project_dir/fives"
barcoded_bams_dir="$project_dir/barcoded_bams_simplex"
methylbed_dir="$project_dir/methylbed"

ref_file=/home/cnorton5/data_hlee308/cnorton5/ref/hg19.fa

module load samtools

#Check for executables
command -v pod5 >/dev/null 2>&1 || { echo "pod5 not found"; exit 1; }
command -v dorado >/dev/null 2>&1 || { echo "dorado not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v modkit >/dev/null 2>&1 || { echo "modkit not found"; exit 1; }

# Check for reference file
[ -f $ref_file ] || { echo "Reference file not found: $ref_file"; exit 1; }

# Make necessary directories
mkdir -p "$project_dir"
mkdir -p "$input_dir"
mkdir -p "$barcoded_bams_dir"
mkdir -p "$methylbed_dir"

shopt -s nullglob
pod5_files=("$input_dir"/*.pod5)
fast5_files=("$input_dir"/*.fast5)

if [ -f "$input_dir/merged.pod5" ]; then
  echo "found merged.pod5 ... continuing with alignment"
elif [ ${#pod5_files[@]} -gt 0 ]; then
  echo "Found pod5 files ... merging pod5 files into merged.pod5"
  pod5 merge "${pod5_files[@]}" --output "$input_dir/merged.pod5" --threads 24
elif [ ${#fast5_files[@]} -gt 0 ]; then
  echo "Found fast5 files ... converting fast5 files to pod5"
  pod5 convert fast5 "${fast5_files[@]}" --output "$input_dir/merged.pod5" --threads 24
else 
  echo "No fast5 or pod5 files found in input directory"
  exit 1
fi

#To get an output name, let's strip the last string after the last "/" in the ref_file
align_genome=$(echo "$ref_file" | awk -F/ '{print $NF}' | awk -F. '{print $1}')
output_name="${align_genome}_simplex.bam"

# Command to run Dorado basecaller with GPUs
dorado basecaller \
  /home/cnorton5/data_hlee308/cnorton5/usr/bin/dorado-0.9.0-linux-x64/bin/models/dna_r9.4.1_e8_hac@v3.3 \
  "$input_dir/merged.pod5" \
  --modified-bases 5mCG_5hmCG \
  --kit-name EXP-NBD104 \
  --reference $ref_file > "$project_dir/$output_name"

# Separate by barcode
dorado demux --output-dir "$barcoded_bams_dir" --no-classify "$project_dir/$output_name"

#rm "$project_dir/$output_name"

# Process each BAM file
for bam_file in "$barcoded_bams_dir"/*.bam; do 

  # Define sorted BAM filename
  sorted_bam_file="${bam_file%.bam}.sorted.bam"

  # Sort the BAM file
  samtools sort --threads 48 -o "$sorted_bam_file" "$bam_file"

  rm "$bam_file"
  
  # Index the sorted BAM file
  samtools index "$sorted_bam_file" -@ 48
  
  # Run modkit on the sorted BAM file
  modkit pileup "$sorted_bam_file" "${bam_file%.bam}.bed" --threads 48
done

#Let's remove the prefix from all files in barcoded_bams_dir, removing everything before the first underscore
for file in "$barcoded_bams_dir"/*; do
  new_name=$(basename "$file" | sed 's/^[^_]*_//')
  new_path="$barcoded_bams_dir/$new_name"

  if [ -e "$new_path" ]; then
    echo "Warning: Skipping rename, '$new_name' already exists."
  else
    mv "$file" "$new_path"
  fi
done

# Move all bed files to methylbed_dir
mv "$barcoded_bams_dir"/*.bed "$methylbed_dir"

echo "Processing complete!"
