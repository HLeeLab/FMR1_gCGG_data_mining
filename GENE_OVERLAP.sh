#!/bin/bash

#SBATCH --partition=shared # 
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=24
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/overlap.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/overlap_err.log # Error log file

set -euo pipefail

# Load required module
module load bedtools

# Define directories
bed_folder="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/hg19_regions_of_interest/hg19_ensembl_annotations"
region_folder="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/hg19_regions_of_interest"
output_folder="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/overlap_hg19/dmr_liu_jaenisch/"

# Create output folder if it doesn't exist
mkdir -p "$output_folder"

# Input files
LIU_BED="$region_folder/CAS9_DMR.bed"
GENES_BED="$bed_folder/hg19_ENSBL_GENES_trimmed.bed"
UP_BED="$bed_folder/hg19_ENSBL_5KB_UP_trimmed.bed"
DOWN_BED="$bed_folder/hg19_ENSBL_5KB_DOWN_trimmed.bed"
EXONS_BED="$bed_folder/hg19_ENSBL_EXONS_trimmed.bed"
INTRONS_BED="$bed_folder/hg19_ENSBL_INTRONS_trimmed.bed"

# Output file
OUTPUT="$output_folder/overlap_matrix.tsv"

# Create an array of unique gene names (using the 4th column)
mapfile -t genes < <(cut -f4 "$GENES_BED" | sort -u)

# Pre-calculate overlaps for each region
tmp_up="$output_folder/tmp_up"
tmp_exons="$output_folder/tmp_exons"
tmp_introns="$output_folder/tmp_introns"
tmp_down="$output_folder/tmp_down"

bedtools intersect -a "$LIU_BED" -b "$UP_BED" -wb | cut -f7 | sort -u > "$tmp_up"
bedtools intersect -a "$LIU_BED" -b "$EXONS_BED" -wb | cut -f7 | sort -u > "$tmp_exons"
bedtools intersect -a "$LIU_BED" -b "$INTRONS_BED" -wb | cut -f7 | sort -u > "$tmp_introns"
bedtools intersect -a "$LIU_BED" -b "$DOWN_BED" -wb | cut -f7 | sort -u > "$tmp_down"

# Print header
echo -e "Gene\tUP\tEXONS\tINTRONS\tDOWN\tANY" > "$OUTPUT"

declare -A up_map exons_map introns_map down_map

# Populate associative arrays with genes found in each region file
while read -r gene; do up_map["$gene"]=1; done < "$tmp_up"
while read -r gene; do exons_map["$gene"]=1; done < "$tmp_exons"
while read -r gene; do introns_map["$gene"]=1; done < "$tmp_introns"
while read -r gene; do down_map["$gene"]=1; done < "$tmp_down"

# Process all genes and check membership in one pass
for gene in "${genes[@]}"; do
    up_val=${up_map["$gene"]:-0}
    exons_val=${exons_map["$gene"]:-0}
    introns_val=${introns_map["$gene"]:-0}
    down_val=${down_map["$gene"]:-0}

    # Calculate "ANY" as 1 if any of the regions has an overlap
    any_val=$(( up_val + exons_val + introns_val + down_val > 0 ? 1 : 0 ))

    echo -e "$gene\t$up_val\t$exons_val\t$introns_val\t$down_val\t$any_val" >> "$OUTPUT"
done

### VERY SLOW GREP
# # Loop over each gene and check membership in each overlap list
# for gene in "${genes[@]}"; do
#     up_val=0
#     exons_val=0
#     introns_val=0
#     down_val=0

#     if grep -q -w "$gene" "$tmp_up"; then
#         up_val=1
#     fi
#     if grep -q -w "$gene" "$tmp_exons"; then
#         exons_val=1
#     fi
#     if grep -q -w "$gene" "$tmp_introns"; then
#         introns_val=1
#     fi
#     if grep -q -w "$gene" "$tmp_down"; then
#         down_val=1
#     fi

#     # Calculate "ANY" as 1 if any of the regions has an overlap
#     any_val=$(( up_val + exons_val + introns_val + down_val > 0 ? 1 : 0 ))

#     echo -e "$gene\t$up_val\t$exons_val\t$introns_val\t$down_val\t$any_val" >> "$OUTPUT"
# done

# Clean up temporary files


rm "$tmp_up" "$tmp_exons" "$tmp_introns" "$tmp_down"

script="/home/cnorton5/data_hlee308/cnorton5/scripts/AddGeneNames_Plot.R"

module load r

##TODO consider adding command line args for the script
Rscript $script

echo "Matrix saved to $OUTPUT"

exit 0



