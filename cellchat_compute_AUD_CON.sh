#!/bin/bash
#SBATCH --job-name=AUD_test                # Job name
#SBATCH --nodes=1                          # Number of nodes
#SBATCH --cpus-per-task=14                 # CPU cores/threads
#SBATCH --mem=800G                         # memory (per node)
#SBATCH --time=0-03:00                     # time (DD-HH:MM)
#SBATCH --partition=zhanglab.p             # use zhanglab partition
#SBATCH --output=output.log                # Standard output
#SBATCH -w galaxy
#SBATCH --mail-type=ALL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jpariser@uci.edu       # Where to send mail


module load conda
source /pkg/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate r8


echo “Running AUD script on galaxy”
Rscript cellchat-AUD_compute.r
Rscript cellchat-CON_compute.r




