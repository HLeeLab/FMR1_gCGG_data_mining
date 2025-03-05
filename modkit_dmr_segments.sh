#!/bin/bash

#SBATCH --partition=shared # 
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=32
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_segment.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_segment_err.log # Error log file

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1

# Set base directories
base_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore
out_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/DMR_Results_hg19

module load samtools

#Check for executables
command -v pod5 >/dev/null 2>&1 || { echo "pod5 not found"; exit 1; }
command -v dorado >/dev/null 2>&1 || { echo "dorado not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v modkit >/dev/null 2>&1 || { echo "modkit not found"; exit 1; }

ref=/home/cnorton5/data_hlee308/cnorton5/ref/hg19.fa

#check that the reference genome is present
[ -f $ref ] || { echo "Reference file not found! $ref"; exit 1; }


dmr_result="$out_dir/single_base_haplotype_dmr.bed"
dmr_segments="$out_dir/single_base_segments.bed"

threads=32

s1="$base_dir/ONT_HG01_hg19/meth_bed/hg19_barcode01.sorted.bed.gz"
s2="$base_dir/ONT_HG01_hg19/meth_bed/hg19_barcode04.sorted.bed.gz"

modkit dmr pair \
  -a ${s1} \
  -b ${s2} \
  -o ${dmr_result} \
  --segment ${dmr_segments} \
  --ref ${ref} \
  --base C \
  --threads ${threads} \
  --force

