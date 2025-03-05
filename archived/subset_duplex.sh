#!/bin/bash

#SBATCH --partition=parallel        # Parallel partition
#SBATCH --time=12:00:00             # Time limit
#SBATCH --nodes=6                # Number of nodes
#SBATCH --ntasks-per-node=24        # Tasks per node (total cores = nodes * ntasks-per-node)
#SBATCH --mem=72G                   # Memory per node
#SBATCH --job-name=pod5_parallel    # Job name
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/subset_output.log  # Standard output
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/subset_error.log    # Standard error
#SBATCH -A hlee308_gpu              # Account for billing

# Load necessary modules
module load samtools
module load parallel

# Verify executables
command -v pod5 >/dev/null 2>&1 || { echo "pod5 not found"; exit 1; }
command -v dorado >/dev/null 2>&1 || { echo "dorado not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "samtools not found"; exit 1; }
command -v modkit >/dev/null 2>&1 || { echo "modkit not found"; exit 1; }
[ -f /home/cnorton5/data_hlee308/cnorton5/ref/T2T-CHM13v2.fna ] || { echo "Reference file not found"; exit 1; }

# Set up directories
base_dir=/home/cnorton5/scr4_hlee308
nanopore_dir="$base_dir/cnorton5/old_nanopore"
pod5_dir="$nanopore_dir/full_pod5"
output_dir="$pod5_dir/split_by_channel"
summary_file="$pod5_dir/summary.tsv"

mkdir -p "$output_dir"


pod5 subset "$pod5_dir" --summary "$pod5_dir/summary.tsv" --columns channel --threads 36 --force-overwrite --output "$pod5_dir/split_by_channel"

# # Extract channels to process
# channels=$(cut -f2 "$summary_file" | tail -n +2 | sort | uniq)

# #Print the first five channels
# echo "First five channels:"
# echo "$channels" | head -n 5


# # Create a task list for parallel processing
# task_list="tasks.txt"
# rm -f "$task_list"
# for channel in $channels; do
#     echo "pod5 subset $pod5_dir --summary $summary_file --columns channel --filter channel $channel --output $output_dir/channel_$channel" >> "$task_list"
# done

# # Use GNU Parallel to distribute tasks
# # Use all available CPUs across all nodes
# total_tasks=$(( SLURM_NTASKS * SLURM_NTASKS_PER_NODE ))
# parallel --jobs "$total_tasks" < "$task_list"

# # Wait for all jobs to complete
# wait

echo "Parallel processing complete."




