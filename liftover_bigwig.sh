#!/bin/bash

#SBATCH --partition=shared # 
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=24
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/liftover_bw.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/liftover_bw.log # Error log file

#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1
base_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/H3K9me3
ref_dir=/home/cnorton5/data_hlee308/cnorton5/ref
INPUT_BW="$base_dir/GSM6754719_1H2-H3K9me3_hg38.bw"
OUTPUT_BW="$base_dir/GSM6754719_1H2-H3K9me3_hg19.bw"

# Input and output file names
CHAIN_FILE="$ref_dir/liftover/hg38ToHg19.over.chain"
HG19_SIZES="$ref_dir/hg19.chrom.sizes"

# Check if input file is provided
if [ -z "$INPUT_BW" ] || [ -z "$OUTPUT_BW" ]; then
    echo "Usage: $0 <input.hg38.bw> <output.hg19.bw>"
    exit 1
fi

# Download required files if not present
if [ ! -f "$CHAIN_FILE" ]; then
    echo "Downloading chain file..."
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
fi

if [ ! -f "$HG19_SIZES" ]; then
    echo "Downloading hg19 chromosome sizes..."
    fetchChromSizes hg19 > hg19.chrom.sizes
fi

# Convert BigWig to BedGraph
echo "Converting BigWig to BedGraph..."
bigWigToBedGraph "$INPUT_BW" $base_dir/temp.hg38.bedGraph

# Run LiftOver
echo "Running LiftOver..."
liftOver $base_dir/temp.hg38.bedGraph "$CHAIN_FILE" $base_dir/hg19.bedGraph unmapped.bed

# Convert BedGraph back to BigWig
echo "Converting BedGraph back to BigWig..."
#First, sort the bedGraph file
sort -k1,1 -k2,2n $base_dir/hg19.bedGraph > $base_dir/hg19.sorted.bedGraph
bedGraphToBigWig $base_dir/hg19.sorted.bedGraph "$HG19_SIZES" "$OUTPUT_BW"

# Cleanup temporary files
echo "Cleaning up..."
rm $base_dir/temp.hg38.bedGraph $base_dir/unmapped.bed $base_dir/hg19.bedGraph

echo "Liftover complete! Output: $OUTPUT_BW"

