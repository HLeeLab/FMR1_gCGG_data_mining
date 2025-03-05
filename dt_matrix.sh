#!/bin/bash

#SBATCH --partition=l40s # Partition with l40s
#SBATCH --time=24:00:00
#SBATCH --nodes=1 # Number of nodes (1 node)
#SBATCH --ntasks-per-node=48
#SBATCH --gres=gpu:4 # 
#SBATCH -A hlee308_gpu # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/dt_matrix.log # Standard output log file
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/dt_matrix_err.log # Error log file


# Set the working directory to the user's home directory 
cd /home/cnorton5 || exit 1


base_dir=/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/ONT_HG01_hg19
methyl_bw="$base_dir/meth_bigwigs"
dt_folder="$base_dir/dt"
regions_of_interest_dir="/home/cnorton5/scr4_hlee308/cnorton5/old_nanopore/ref/hg19_regions_of_interest"
regions_name="Liu_Jaenisch_FXS_supplemental1"
regions_path="$regions_of_interest_dir/$regions_name.bed"
mod_type="hmc"


mkdir -p $dt_folder

module load anaconda 
conda activate dt_env

computeMatrix scale-regions -S "$methyl_bw"/*."$mod_type".bw -R "$regions_path" -o "$dt_folder/$regions_name.$mod_type.matrix.gz"
plotHeatmap -m "$dt_folder/$regions_name.$mod_type.matrix.gz" -out "$dt_folder/$regions_name.$mod_type.heatmap.png" --colorMap RdBu --plotTitle "Heatmap of methylation levels in $regions_name" --plotFileFormat png
conda deactivate






