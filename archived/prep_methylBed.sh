#!/bin/bash

#SBATCH --partition=shared # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=24 # Number of tasks per node (24 cores)
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/prep_methylBed_output.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/prep_methylBed_error.log # Error log file

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1

# Set base directories
base_dir=/home/cnorton5/scr4_hlee308
nanopore_dir="$base_dir/cnorton5/old_nanopore/ONT_HG02_hg19_full"
methylbed_dir="$nanopore_dir/methylbed"


module load samtools
module load tabix
module load bedtools

#Check for executables
command -v pod5 >/dev/null 2>&1 || { echo "pod5 not found"; exit 1; }
command -v dorado >/dev/null 2>&1 || { echo "dorado not found"; exit 1; }
command -v modkit >/dev/null 2>&1 || { echo "modkit not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "bgzip not found"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "tabix not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found"; exit 1; }

#First, unzip all the .bed.gz files in methyl_bed
for file in "$methylbed_dir"/*.bed.gz; do
    gunzip "$file"
done

#Let's zip and index the .bed files in methyl_bed
for file in "$methylbed_dir"/*.bed; do
    bgzip -c "$file" > "$file".gz -@ 24
    tabix -p bed "$file".gz
done

