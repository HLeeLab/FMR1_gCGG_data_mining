#!/bin/bash

#SBATCH --partition=shared # 
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=24
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/temp.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/temp.log # Error log file

#!/bin/bash
# Change to the working directory
cd /home/cnorton5 || exit 1

module load bedtools

# Set output directory
OUTDIR="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/ref/hg19_regions_of_interest/hg19_ensembl_annotations"
mkdir -p "$OUTDIR"

# Download hg19 ENSEMBL gene annotations from UCSC Table Browser
echo "Downloading ENSEMBL gene annotations for hg19..."
wget -O "$OUTDIR/hg19_ensembl_genes.bed" "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz"

# Unzip file
gunzip -c "$OUTDIR/hg19_ensembl_genes.bed" > "$OUTDIR/hg19_ensembl_genes.txt"

# Extract columns: chr, txStart, txEnd, strand, exon starts/ends, gene name
awk 'BEGIN {OFS="\t"} {print $3, $5, $6, $2, $4, $10, $11}' "$OUTDIR/hg19_ensembl_genes.txt" > "$OUTDIR/genes.bed"

# Extract 5' UTR (5KB upstream) and 3' UTR (5KB downstream)
echo "Extracting 5' UTR and 3' UTR..."
awk -v out="$OUTDIR" 'BEGIN {OFS="\t"} {
    if ($5 == "+") {
        print $1, ($2 - 5000 < 0 ? 0 : $2 - 5000), $2, $4, $5 > out"/5p_UTR_5kb.bed";
        print $1, $3, $3 + 5000, $4, $5 > out"/3p_UTR_5kb.bed";
    } else {
        print $1, $3, $3 + 5000, $4, $5 > out"/5p_UTR_5kb.bed";
        print $1, ($2 - 5000 < 0 ? 0 : $2 - 5000), $2, $4, $5 > out"/3p_UTR_5kb.bed";
    }
}' "$OUTDIR/genes.bed"

# Extract exons
echo "Extracting exons..."
awk 'BEGIN {OFS="\t"} {
    split($6, starts, ",");
    split($7, ends, ",");
    for (i=1; i<=length(starts); i++) {
        if (starts[i] != "" && ends[i] != "") {
            print $1, starts[i], ends[i], $4, $5
        }
    }
}' "$OUTDIR/genes.bed" | sort -k1,1 -k2,2n | uniq > "$OUTDIR/exons.bed"

# Extract introns using bedtools
echo "Extracting introns..."
bedtools subtract -a <(sort -k1,1 -k2,2n "$OUTDIR/genes.bed") -b <(sort -k1,1 -k2,2n "$OUTDIR/exons.bed") > "$OUTDIR/introns.bed"



echo "Processing complete! Results are saved in $OUTDIR"
