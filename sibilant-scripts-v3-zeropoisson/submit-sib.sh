#!/bin/bash -l
# SJ Riesenfeld & J Farrell

## TO RUN THIS SCRIPT:
## sbatch submit-sib.sh [directory_name]

## ===================================
## SET BSUB JOB OPTIONS APPROPRIATELY!
## ===================================

#SBATCH -J sibilant
#SBATCH -n 16              ## num threads;
#SBATCH -N 1               ## All cores on one node.
#SBATCH --mem-per-cpu=3072        ## 4G of RAM per node
#SBATCH -t 0-04:00         ## Time job can run (4 hours)
#SBATCH -p serial_requeue  ## Queue
#SBATCH -o logs/UMI_PL.%A.out
#SBATCH -e logs/UMI_PL.%A.err

# ============================
## SET JOB-SPECIFIC VARIABLES!
# ============================

function v_exe
{
    echo "$1"
    eval "$1" || error_exit "Cannot execute command: $1"
}

# Load package
source new-modules.sh
module load matlab/R2016b-fasrc01
module load java/1.8.0_45-fasrc01 

matlab -nodisplay -r "sample_dir = '${1}'; run('start_sibilant.m'); exit;"