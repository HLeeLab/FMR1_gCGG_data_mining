#!/bin/bash

#SBATCH --partition=l40s # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=48 # Number of tasks per node (24 cores)
#SBATCH --mem=96G # Memory allocation (96 GB)
#SBATCH --gres=gpu:4 # 
#SBATCH -A hlee308_gpu # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_dmr.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_dmr_err.log # Error log file

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1

# Set base directories
base_dir=/home/cnorton5/scr4_hlee308
nanopore_dir="$base_dir/cnorton5/old_nanopore"
barcoded_bams_dir="$nanopore_dir/barcoded_bams_simplex"
methylbed_dir="$barcoded_bams_dir/methylbed"
results_dir="$barcoded_bams_dir/dmrResults"

module load samtools

#Check for executables
command -v pod5 >/dev/null 2>&1 || { echo "pod5 not found"; exit 1; }
command -v dorado >/dev/null 2>&1 || { echo "dorado not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v modkit >/dev/null 2>&1 || { echo "modkit not found"; exit 1; }


#First, check that there are cpg islands
[ -f /home/cnorton5/data_hlee308/cnorton5/ref/T2T_other/T2T_cpg_island_refseq_simple.bed ] || { echo "Reference file not found: $base_dir/ref//T2T_other/T2T_cpg_island_refseq_simple.bed"; exit 1; }

#check that the reference genome is present
[ -f /home/cnorton5/data_hlee308/cnorton5/ref/T2T-CHM13v2.fna ] || { echo "Reference file not found: $base_dir/ref/T2T_other/T2T.fa"; exit 1; }

# Make necessary directories
mkdir -p "$results_dir"

regions=/home/cnorton5/data_hlee308/cnorton5/ref/T2T_other/T2T_cpg_island_refseq_simple.bed
dmr_result="$results_dir/output.bed"
ref=/home/cnorton5/data_hlee308/cnorton5/ref/T2T-CHM13v2.fna
threads=48

modkit dmr pair --ref ${ref} -a "${methylbed_dir}/EXP-NBD104_barcode01.bed.gz" -b "${methylbed_dir}/EXP-NBD104_barcode04.bed.gz" -o ${dmr_result} -r ${regions} --base C --threads ${threads}
