#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -A hlee308
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/ChIP_seq_output.log
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/ChIP_seq_error.log

# Change to the working directory
cd /home/cnorton5 || exit 1

module load samtools
module load fastqc
module load bowtie2
module load anaconda


conda activate chipseq
# Set base directory for results
base_dir="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/ChIP_Seq_hg19"
ref_dir="/home/cnorton5/data_hlee308/cnorton5/ref"

mkdir -p "$base_dir/bowtie2_results"
mkdir -p "$base_dir/qc_results"
mkdir -p "$base_dir/trimmed_reads"

BOWTIE2_INDEX="$ref_dir/hg19_other/hg19_bowtie2_index/hg19_bowtie2_index"

# Directory containing raw FASTQ files
IN_DIR="/home/cnorton5/data_hlee308/erisOne/fastq_hg19"
TRIM_DIR="$base_dir/trimmed_reads"
OUT_DIR="$base_dir/bowtie2_results"

mkdir -p "$OUT_DIR"

# Define sample names
declare -a SAMPLES=(
    "hFXS-iPS_Cas9-ChIP-line1.GSM2742477"
    "hFXS-iPS_Cas9-ChIP-line2.GSM2742480"
    "hFXS-iPS_Cas9-ChIP-line3.GSM2742481"
)

mock="hFXS-iPS_Cas9-ChIP-MockLine.GSM2742479"

# Add mock sample
SAMPLES+=("$mock")

# Preprocessing and alignment
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing $SAMPLE..."

    R1="${IN_DIR}/${SAMPLE}.fq.gz"
    TRIMMED_R1="${TRIM_DIR}/${SAMPLE}_trimmed.fq.gz"


    trim_galore --gzip --cores 8 --fastqc -o "$TRIM_DIR" "$R1"

    # Create output directory for sample
    SAMPLE_OUT="$OUT_DIR/$SAMPLE"
    mkdir -p "$SAMPLE_OUT"

    # Align with Bowtie2
    bowtie2 -p 12 -x "$BOWTIE2_INDEX" -U "$TRIMMED_R1" 2> "$SAMPLE_OUT/${SAMPLE}_bowtie2.log" | \
        samtools view -bS - | \
        samtools sort -@ 12 -o "$SAMPLE_OUT/${SAMPLE}_sorted.bam"

done

#Remove mock from SAMPLES array
SAMPLES=("${SAMPLES[@]/$mock}")


#Call peaks for each sample, except the mock, using MACS3
module load MACS
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_OUT="$OUT_DIR/$SAMPLE"
    MOCK_OUT="$OUT_DIR/$mock"

    echo "Calling peaks for $SAMPLE..."
    macs3 callpeak -t "$SAMPLE_OUT/${SAMPLE}_sorted.bam" -c "$MOCK_OUT/${mock}_sorted.bam" -f BAM -g hs -n "$SAMPLE" --broad --outdir "$SAMPLE_OUT"
done
module unload MACS

conda deactivate

