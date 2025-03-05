#!/bin/bash

#SBATCH --partition=shared # 
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=24
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/temp.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/temp.log # Error log file

module load samtools
module load htslib



base_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/hg19/barcoded_bams_simplex


for file in "$base_dir"/*.bed; do
    bgzip -@ 12 -c "$file" > "${file}.gz"  # Create compressed file
    tabix -p bed "${file}.gz"        # Index the compressed file
done