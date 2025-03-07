#!/bin/bash
#SBATCH --partition=shared        # Use the l40s partition (CPU nodes)
#SBATCH --time=24:00:00             # Job time limit
#SBATCH --nodes=1             # Run on a single node
#SBATCH --ntasks-per-node=24     # Use 24 CPU cores
#SBATCH -A hlee308            
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/hisat2_output.log  # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/hisat2_error.log    # Error log file

# Change to the working directory
cd /home/cnorton5 || exit 1

module load samtools
module load fastqc
module load anaconda

# Set base directory for results (modify as needed)
base_dir="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/RNA_HG01_hg19"
ref_dir="/home/cnorton5/data_hlee308/cnorton5/ref"

HISAT2_INDEX="/home/cnorton5/data_hlee308/cnorton5/ref/hg19_other/hg19_hisat2_index/hg19"

# Directory containing FASTQ files
IN_DIR="$base_dir/HG02mRNA"
# Output directory for HISAT2 alignments
OUT_DIR="$base_dir/hisat2_hg19"
mkdir -p "$OUT_DIR"

# Define sample names (modify these names to match your FASTQ file naming)
declare -a SAMPLES=(
    "HG02mRNA-FXS_gNHG3_36d_1"
    "HG02mRNA-FXS_gNHG3_36d_2"
)

#This is the conda environment that has trim_galore installed
conda activate chipseq


# Loop over each sample and run HISAT2
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing $SAMPLE..."

    # Define input FASTQ file paths (assumes paired-end reads with _R1 and _R2 suffixes)
    R1="${IN_DIR}/${SAMPLE}_R1.fq.gz"
    R2="${IN_DIR}/${SAMPLE}_R2.fq.gz"

    trimmed_R1="${IN_DIR}/${SAMPLE}_R1_val_1.fq.gz"
    trimmed_R2="${IN_DIR}/${SAMPLE}_R2_val_2.fq.gz"

    trim_galore --paired --gzip --fastqc --cores 8 -o "$IN_DIR" "$R1" "$R2"

    # Create an output subdirectory for the sample
    SAMPLE_OUT="$OUT_DIR/$SAMPLE"
    mkdir -p "$SAMPLE_OUT"
    
    # Option 2 (Recommended): Pipe directly to SAMtools to produce a sorted BAM file.
    
    hisat2 -p 12 -x "$HISAT2_INDEX" -1 "$trimmed_R1" -2 "$trimmed_R2" 2> "$SAMPLE_OUT/${SAMPLE}_hisat2.log" | \
        samtools view -bS - | \
        samtools sort -@ 12 -o "$SAMPLE_OUT/${SAMPLE}_sorted.bam"
    
    # Then you can index the BAM file:
    samtools index "$SAMPLE_OUT/${SAMPLE}_sorted.bam"

done

conda deactivate

