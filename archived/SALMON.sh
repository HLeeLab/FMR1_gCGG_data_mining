#!/bin/bash

#SBATCH --partition=l40s # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=24 # Number of tasks per node (24 cores)
#SBATCH --mem=48G # Memory allocation (96 GB)
#SBATCH --gres=gpu:2 # 
#SBATCH -A hlee308_gpu # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/salmon_output.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/salmon_error.log # Error log file

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1

# RNA-seq preprocessing: Quantification with Salmon
# Ensure Salmon is installed and available in PATH

base_dir="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/hg19/rna_seq"

mkdir -p $base_dir/salmon_results

# Define paths
SALMON_INDEX="/home/cnorton5/data_hlee308/cnorton5/ref/hg19_other/hg19_salmon_k31"  # Path to prebuilt Salmon index
IN_DIR="/home/cnorton5/data_hlee308/erisOne/data_HG-NGS/HG01mRNA" # Directory containing FASTQ files
OUT_DIR="$base_dir/salmon_results" # Output directory for quantification

# Create output directory
mkdir -p $OUT_DIR

# List of sample FASTQ pairs
declare -a SAMPLES=(
    "HG01mRNA-CGGonly_24d_mT_1"
    "HG01mRNA-CGGonly_24d_mT_2"
    "HG01mRNA-dCas9E_CGG_1"
    "HG01mRNA-dCas9E_CGG_2"
)

# Quantification loop
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing $SAMPLE..."

    # Define input FASTQ files
    R1="${IN_DIR}/${SAMPLE}_R1.fq.gz"
    R2="${IN_DIR}/${SAMPLE}_R2.fq.gz"
    SAMPLE_OUT="$OUT_DIR/$SAMPLE"
    SAMPLE_OUT_SAM="$SAMPLE_OUT/$SAMPLE.sam"

    
    # Run Salmon quantification
    salmon quant -l A -i $SALMON_INDEX -1 $R1 -2 $R2 --validateMappings --gcBias --threads 12 -o "$SAMPLE_OUT"  #
    #--writeMappings="$SAMPLE_OUT_SAM"

done

