#!/bin/bash

#SBATCH --partition=shared # 
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=24
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/temp.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/temp.log # Error log file

#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Input and output file names
INPUT_BW="$1"
OUTPUT_BW="$2"
CHAIN_FILE="hg38ToHg19.over.chain.gz"
HG19_SIZES="hg19.chrom.sizes"

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
bigWigToBedGraph "$INPUT_BW" temp.hg38.bedGraph

# Run LiftOver
echo "Running LiftOver..."
liftOver temp.hg38.bedGraph "$CHAIN_FILE" temp.hg19.bedGraph unmapped.bed

# Convert BedGraph back to BigWig
echo "Converting BedGraph back to BigWig..."
bedGraphToBigWig temp.hg19.bedGraph "$HG19_SIZES" "$OUTPUT_BW"

# Cleanup temporary files
echo "Cleaning up..."
rm temp.hg38.bedGraph temp.hg19.bedGraph unmapped.bed

echo "Liftover complete! Output: $OUTPUT_BW"

