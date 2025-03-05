#!/bin/bash
#SBATCH --partition=shared
#SBATCH --time=02:00:00  # Set to 2 hours
#SBATCH --nodes=1        # Use a single node
#SBATCH --ntasks-per-node=32 # Request 64 cores
#SBATCH -A hlee308       # Account for billing
#SBATCH --output=/home/cnorton5/data_hlee308/cnorton5/logs/liftover_output.log  # Standard output log
#SBATCH --error=/home/cnorton5/data_hlee308/cnorton5/logs/liftover_error.log  # Error log

# Load R module (modify if necessary for your cluster)
module load r

# Run the liftover script
Rscript /home/cnorton5/data_hlee308/cnorton5/scripts/liftover.R

