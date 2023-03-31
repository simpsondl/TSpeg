#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --job-name="makeinters"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out

# Allows batch submission of prepare-intermediate_files.sh

# Call script with                 
### sbatch --export=indir=INDIR prepare_intermediate_files_driver.s
### INDIR is full path for input directory containing bam files

module purge

/Genomics/grid/users/ds65/scripts/pegRNA/prepare_intermediate_files.sh $indir
