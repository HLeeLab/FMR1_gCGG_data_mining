#!/bin/bash

#SBATCH --partition=shared # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=24
#SBATCH -A hlee308_gpu # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_region_stats.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_region_stats_err.log # Error log file


base_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/
region_stats_dir="$base_dir/DMR_Results_hg19"

regions_of_interest_dir="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/hg19_regions_of_interest"

region_file_names=("Liu_Jaenisch_FXS_supplemental1.bed")


region_file_paths=()

#Sometimes, there are name, score, and strand columns in the region files. We want to keep these. 
#TODO add a check to verify that columns 4, 5, and 6 are name, score, and strand columns
for region_file_name in "${region_file_names[@]}"; do
    region_file="$regions_of_interest_dir/$region_file_name"
    if [[ -f "$region_file" ]]; then
        region_file_paths+=("$region_file")
    fi
done

# Make necessary directories
mkdir -p "$region_stats_dir"

#Let's declare a list of all the methyl bed files we want to process
methyl_bed_file_paths=()
methyl_bed_file_paths+=("$base_dir/ONT_HG01_hg19/methylbed/hg19_barcode01.sorted.bed.gz")
methyl_bed_file_paths+=("$base_dir/ONT_HG01_hg19/methylbed/hg19_barcode04.sorted.bed.gz")
methyl_bed_file_paths+=("$base_dir/ONT_HG02_hg19_full/methylbed/EXP-NBD104_barcode01.bed.gz")
methyl_bed_file_paths+=("$base_dir/ONT_HG02_hg19_full/methylbed/EXP-NBD104_barcode04.bed.gz")


for input_file in "${methyl_bed_file_paths[@]}"; do
    echo "Processing $input_file..."

    for region_file in "${region_file_paths[@]}"; do
        region_file_name=$(basename "$region_file")
        file_name=$(basename "$input_file")
        output_file="${region_stats_dir}/${file_name}_${region_file_name}.bed_stats.tsv"
        modkit stats ${input_file} --regions ${region_file} -o "$output_file" --mod-codes "h,m" --threads 24
    done

    echo "Finished processing barcode ${input_file}."
done
