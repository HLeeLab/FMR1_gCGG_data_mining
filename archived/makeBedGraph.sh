#!/bin/bash

# Loop through barcode files

cd /Users/carternorton/Desktop/barcoded_bams_simplex
mkdir -p bedgraphs
for barcode in 04; do
    # Define the input file name
    input_file="EXP-NBD104_barcode${barcode}.bed"

    # Check if the file exists
    if [[ -f "$input_file" ]]; then
        echo "Processing $input_file..."

        # Extract lines with 'h' in the 4th column into a temporary file
        awk '$4 == "h"' "$input_file" > temp_h.bed

        # Extract lines with 'm' in the 4th column into a temporary file
        awk '$4 == "m"' "$input_file" > temp_m.bed

        # Process the 'h' file and save the output
        awk -v OFS="\t" '{print $1, $2, $3, $11, $10}' temp_h.bed > "bedgraphs/EXP-NBD104_barcode${barcode}_h.bed"

        # Process the 'm' file and save the output
        awk -v OFS="\t" '{print $1, $2, $3, $11, $10}' temp_m.bed > "bedgraphs/EXP-NBD104_barcode${barcode}_m.bed"

        # Remove temporary files
        rm temp_h.bed temp_m.bed

        echo "Finished processing barcode ${barcode}."
    else
        echo "File $input_file does not exist, skipping..."
    fi
done
