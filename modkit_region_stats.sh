#!/bin/bash

#SBATCH --partition=l40s # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=48
#SBATCH --gres=gpu:4 # 
#SBATCH -A hlee308_gpu # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_region_stats.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/modkit_region_stats_err.log # Error log file


base_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/ONT_HG01_hg19
methyl_bed_dir="$base_dir/meth_bed"
region_stats_dir="$base_dir/region_stats"

regions_of_interest_dir="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/ref/hg19_regions_of_interest"

region_file_names=("hg19_GENCODE_1000bp_TSS.bed")


region_file_paths=()

#Use a for loop to create region file paths and if they have more than 3 columns, trim them down to 3 and rename them as "shortened_" + region_file_name
# for region_file_name in "${region_file_names[@]}"; do
#     region_file="$regions_of_interest_dir/$region_file_name"
#     if [[ -f "$region_file" ]]; then
#         if [[ $(head -n 1 "$region_file" | awk '{print NF}') -gt 3 ]]; then
#             awk -v OFS='\t' '{print $1, $2, $3}' "$region_file" > "$regions_of_interest_dir/shortened_$region_file_name"
#             region_file_paths+=("$regions_of_interest_dir/shortened_$region_file_name")
#         else
#             region_file_paths+=("$region_file")
#         fi
#     fi
# done

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

# First, check if there are any files ending in .bed.gz in the methyl_bed_dir, if not, check for .bed files
if [[ $(ls $methyl_bed_dir | grep ".bed.gz" | wc -l) -eq 0 ]]; then
    methyl_bed_files=$(ls $methyl_bed_dir | grep ".bed")
    if [[ $(ls $methyl_bed_dir | grep ".bed" | wc -l) -eq 0 ]]; then
        echo "No .bed or .bed.gz files found in $methyl_bed_dir, exiting..."
        exit 1
    else
        for input_file in $methyl_bed_files; do
            echo "Compressing $input_file... and indexing"
            bgzip -c "$methyl_bed_dir/$input_file" > "$methyl_bed_dir/$input_file.gz"
            tabix -p bed "$methyl_bed_dir/$input_file.gz"
        done
    fi
else
    echo "Found .bed.gz files in $methyl_bed_dir."
fi

#Let's get a list of all files in the methyl_bed_dir that end in .bed.gz
methyl_bed_file_paths=("$methyl_bed_dir"/*.bed.gz)

#Print all the files in the methyl_bed_dir
echo "Files in $methyl_bed_dir:"
for file in "${methyl_bed_file_paths[@]}"; do
    echo "$file"
done



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
