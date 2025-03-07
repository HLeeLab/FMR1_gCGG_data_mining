#!/bin/bash

#SBATCH --partition=shared # 
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=32
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_multi.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_multi_err.log # Error log file

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1

# Set base directories
base_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore
out_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/DMR_Results_hg19/multi

mkdir -p $out_dir

module load samtools

#Check for executables
command -v pod5 >/dev/null 2>&1 || { echo "pod5 not found"; exit 1; }
command -v dorado >/dev/null 2>&1 || { echo "dorado not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v modkit >/dev/null 2>&1 || { echo "modkit not found"; exit 1; }

ref=/home/cnorton5/data_hlee308/cnorton5/ref/hg19.fa

#check that the reference genome is present
[ -f $ref ] || { echo "Reference file not found! $ref"; exit 1; }


regions="$base_dir/hg19_regions_of_interest/Liu_Jaenisch_FXS_supplemental1.bed"
dmr_result="$out_dir/Liu_Jaenisch_dmr.bed"
ind_seg_result="$out_dir/1v4_ind_segment.bed"
seg_result="$out_dir/1v4_segment.bed"

threads=32

s1="$base_dir/ONT_HG01_hg19/methylbed/hg19_barcode01.sorted.bed.gz"
s1a="$base_dir/ONT_HG02_hg19_full/methylbed/EXP-NBD104_barcode01.bed.gz"
s2="$base_dir/ONT_HG01_hg19/methylbed/hg19_barcode04.sorted.bed.gz"
s2a="$base_dir/ONT_HG02_hg19_full/methylbed/EXP-NBD104_barcode04.bed.gz"



modkit dmr multi \
  -s ${s1} cont1 \
  -s ${s1a} cont2\
  -s ${s2} treat1\
  -s ${s2a} treat2\
  -o ${dmr_result} \
  --regions ${regions} \
  --ref ${ref} \
  --base C \
  --threads ${threads} \
  --force


modkit dmr pair \
  -a ${s1} \
  -a ${s1a} \
  -b ${s2} \
  -b ${s2a} \
  -o ${ind_seg_result} \
  --segment ${seg_result} \
  --ref ${ref} \
  --base C \
  --threads ${threads} \
  --force

