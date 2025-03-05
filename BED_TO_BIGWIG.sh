#!/bin/bash

#SBATCH --partition=l40s # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=48 # Number of tasks per node (24 cores)
#SBATCH --mem=96G # Memory allocation (96 GB)
#SBATCH --gres=gpu:4 # 
#SBATCH -A hlee308_gpu # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/BED_TO_BIGWIG.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/BED_TO_BIGWIG_error.log # Error log file

# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1

# Set base directories
base_dir=/home/cnorton5/scr4_hlee308/
input_dir="$base_dir/cnorton5/old_nanopore/ONT_HG01_hg19/barcoded_bams_simplex"
project_dir="$base_dir/cnorton5/old_nanopore/ONT_HG01_hg19"
output_dir="$base_dir/cnorton5/old_nanopore/ONT_HG01_hg19/meth_bigwigs"

chrom_sizes="/home/cnorton5/data_hlee308/cnorton5/ref/hg19.chrom.sizes"

mkdir -p $output_dir

#Identify all bed files in the input directory
bed_files=$(find $input_dir -name "*.bed")

#Iterate over all bed files
for bed_file in $bed_files; do

    #Let's create a temporory bed_file that only includes column1 values that are located in the first column of the chrom_sizes file
    awk 'NR==FNR{a[$1];next}($1 in a)' $chrom_sizes $bed_file > "$output_dir/temp.bed"

    #Sort the bed file by chromosome
    sort -k1,1 -k2,2n "$output_dir/temp.bed" > "$output_dir/temp.sorted.bed"

    #Overwrite the temp.bed file with the sorted bed file
    mv "$output_dir/temp.sorted.bed" "$output_dir/temp.bed"

    echo "Processing $bed_file"
    #Extract the sample name from the bed file
    sample_name=$(basename $bed_file .bed)

    #convert into two temp files, depending on hmc or mC 
    #modkit bedmethyl tobigwig  --sizes $chrom_sizes --mod-codes h --nthreads 12 $bed_file "$output_dir/$sample_name.hmc.bw"
    #modkit bedmethyl tobigwig  --sizes $chrom_sizes --mod-codes m --nthreads 12 $bed_file "$output_dir/$sample_name.mc.bw"

    modkit bedmethyl tobigwig  --sizes $chrom_sizes --mod-codes h --nthreads 12 "$output_dir/temp.bed" "$output_dir/$sample_name.hmc.bw"
    modkit bedmethyl tobigwig  --sizes $chrom_sizes --mod-codes m --nthreads 12 "$output_dir/temp.bed" "$output_dir/$sample_name.mc.bw"

    #Remove the temp file
    rm "$output_dir/temp.bed"


    # #Now we are going to try to a more roundabout method 
    # #Separate into two files, one for hmc and one for mc
    # awk '{if($4 == "h") print $0}' $bed_file > "$output_dir2/$sample_name.hmc.bed"
    # awk '{if($4 == "m") print $0}' $bed_file > "$output_dir2/$sample_name.mc.bed"

    # #Now convert to bedgraph by keeping columns 1,2,3,10
    # cut -f1,2,3,10 "$output_dir2/$sample_name.hmc.bed" > "$output_dir2/$sample_name.hmc.bedgraph"
    # cut -f1,2,3,10 "$output_dir2/$sample_name.mc.bed" > "$output_dir2/$sample_name.mc.bedgraph"

    # #Now convert to bigwig using bedGraphToBigWig
    # bedGraphToBigWig "$output_dir2/$sample_name.hmc.bedgraph" $chrom_sizes "$output_dir2/$sample_name.hmc.bw"
    # bedGraphToBigWig "$output_dir2/$sample_name.mc.bedgraph" $chrom_sizes "$output_dir2/$sample_name.mc.bw"
    
    #Now remove the temp files
    # rm "$output_dir2/$sample_name.hmc.bed"
    # rm "$output_dir2/$sample_name.mc.bed"
    # rm "$output_dir2/$sample_name.hmc.bedgraph"
    # rm "$output_dir2/$sample_name.mc.bedgraph"
done





