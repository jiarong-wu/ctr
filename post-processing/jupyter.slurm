#!/bin/bash
#SBATCH --job-name="Jupyter"                # Job name
#SBATCH --mail-user=jiarongw@stanford.edu  # Email address    
#SBATCH --mail-type=NONE                    # Mail notification type (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --partition=cpu                     # Node partition
#SBATCH --nodes=1                           # Number of nodes requested
#SBATCH --ntasks=8                          # Number of processes
#SBATCH --time=14:00:00                     # Time limit request

date
conda init
conda activate ww
EXEC_DIR=$HOME
hostname && jupyter-notebook --no-browser --notebook-dir=$EXEC_DIR
