#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --job-name="makeinters"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out

##########################################################################
#### Before running this, edit the SBATCH directives above to point   ####
#### to your own log directory and to have appropriate time and       #### 
#### memory constraints based on the amount of sequencing performed   #### 
#### and number of sequences that need to be aligned.                 ####
##########################################################################

##########################################################################
#### This script expects a directory to be given as input. Two sorted ####
#### BAM files, resulting from running bowtie2_single_read_align.s on ####
#### input FASTQs (forward and reverse), should be inside the input   #### 
#### directory. This script will then check reads for recombination-- ####
#### where the pegRNA does not match the target site--and generate a  #### 
#### whitelist of reads to move forward in addition to coverage stats ####
#### for all pegRNAs in the library.                                  ####
##########################################################################

##########################################################################
#### Call script with                                                 ####
#### sbatch --export=indir=INDIR prepare_intermediate_files_driver.s  ####
#### where                                                            ####
#### #### INDIR is full path for input directory containing bam files ####
#### #### #### The basename of INDIR will appear in file outputs.     ####
##########################################################################

module purge

/Genomics/grid/users/ds65/scripts/pegRNA/prepare_intermediate_files.sh $indir
