#!/bin/bash

#SBATCH --partition=shared # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=32
#SBATCH -A hlee308 # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/dt_matrix.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/dt_matrix_err.log # Error log file


# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1


base_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore
dt_folder="$base_dir/dt"
regions_of_interest_dir="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/hg19_regions_of_interest"
regions_name="Liu_Jaenisch_FXS_supplemental1"
regions_path="$regions_of_interest_dir/$regions_name.bed"
mod_type="mc"


mkdir -p $dt_folder

module load anaconda 
conda activate dt_env

#Let's grab each of the bam files of interest, by declaring a list of files

barcode1a="$base_dir/ONT_HG01_hg19/meth_bigwigs/hg19_barcode01.sorted.mc.bw"
barcode1b="$base_dir/ONT_HG02_hg19/meth_bigwigs/EXP-NBD104_barcode01.mc.bw"
barcode4a="$base_dir/ONT_HG01_hg19/meth_bigwigs/hg19_barcode04.sorted.mc.bw"
barcode4b="$base_dir/ONT_HG02_hg19/meth_bigwigs/EXP-NBD104_barcode04.mc.bw"

#Let's make a list to pass to computeMatrix
bw_list="$barcode1a $barcode1b $barcode4a $barcode4b"

computeMatrix scale-regions -S $bw_list -R "$regions_path" -o "$dt_folder/$regions_name.$mod_type.matrix.gz" -p 24
plotHeatmap -m "$dt_folder/$regions_name.$mod_type.matrix.gz" -out "$dt_folder/$regions_name.$mod_type.heatmap.png" --colorMap RdBu --plotTitle "Heatmap of methylation levels in $regions_name" --plotFileFormat png --kmeans 4
conda deactivate






