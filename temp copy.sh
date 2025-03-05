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


#!/bin/bash

# Input files
LIU_BED="Liu_Jaenisch_FXS_supplemental1.bed"
GENES_BED="hg19_ENSBL_GENES_trimmed.bed"
UP_BED="hg19_ENSBL_5KB_UP_trimmed.bed"
DOWN_BED="hg19_ENSBL_5KB_DOWN_trimmed.bed"
EXONS_BED="hg19_ENSBL_EXONS_trimmed.bed"
INTRONS_BED="hg19_ENSBL_INTRONS_trimmed.bed"

# Output file
OUTPUT="overlap_matrix.tsv"

# Extract gene names from ENSBL_GENES
genes=$(cut -f4 "$GENES_BED" | sort | uniq)

# Print header
echo -e "Gene\tUP\tEXONS\tINTRONS\tDOWN\tANY" > "$OUTPUT"

# Function to check for overlap
check_overlap() {
    local bed_file=$1
    local gene_list=$2
    local tmp_file=$(mktemp)

    # Find overlapping genes
    bedtools intersect -a "$GENES_BED" -b "$bed_file" | cut -f4 | sort | uniq > "$tmp_file"

    # Create associative array for quick lookup
    declare -A overlap
    while read -r gene; do
        overlap["$gene"]=1
    done < "$tmp_file"

    # Print result for each gene
    for gene in $gene_list; do
        if [[ -n "${overlap[$gene]}" ]]; then
            echo 1
        else
            echo 0
        fi
    done

    rm "$tmp_file"
}

# Compute overlap for each region
up=$(check_overlap "$UP_BED" "$genes")
exons=$(check_overlap "$EXONS_BED" "$genes")
introns=$(check_overlap "$INTRONS_BED" "$genes")
down=$(check_overlap "$DOWN_BED" "$genes")

# Construct the matrix
paste <(echo "$genes") <(echo "$up") <(echo "$exons") <(echo "$introns") <(echo "$down") | while IFS=$'\t' read -r gene u e i d; do
    any=$((u + e + i + d > 0 ? 1 : 0))
    echo -e "$gene\t$u\t$e\t$i\t$d\t$any"
done >> "$OUTPUT"

echo "Matrix saved to $OUTPUT"
