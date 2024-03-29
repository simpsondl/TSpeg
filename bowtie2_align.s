#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --job-name="align"
#SBATCH --output=%x-%j.out

##########################################################################
#### BEFORE RUNNING THIS, edit the SBATCH directives above to point   ####
#### to your own log directory and to have appropriate time and       ####
#### memory constraints based on the amount of sequencing performed   ####
#### and number of sequences that need to be aligned. Additionally,   ####
#### check the DEPENDENCY section below for scripts that are called   ####
#### by this file, then adjust file paths to point to a local copy    ####
#### and make any necessary changes to those scripts themselves.      ####
##########################################################################

##########################################################################
#### PURPOSE: Align each read independently to reference              ####
##########################################################################

##########################################################################
#### PREREQUISITE: NONE                                               ####
##########################################################################

##########################################################################
#### DEPENDENCY: NONE                                                 ####
##########################################################################

##########################################################################
#### OUTPUT: Creates a sorted BAM file                                ####
##########################################################################

##########################################################################
#### This script expects a single fastq file, either R1 or R2, and    ####
#### will align that fastq file to the indicated reference, provided  ####
#### by the index variable in the submission call. After aligning all #### 
#### reads, the resulting alignment is sorted according to read names ####
#### and saved in BAM format.                                         ####
##########################################################################

###############################################################################################
#### Call script with                                                                      ####
#### sbatch --export=in=READ.FQ,index=REF,out=OUT bowtie2_align.s                          ####
#### where                                                                                 ####
#### #### READ.FQ is the input fastq file to be processed (full path, can be gz)           ####
#### #### REF is the base name of the reference index (full path)                          ####
#### #### OUT is the name of the output file (full path, should include .bam as extension) ####
###############################################################################################

module purge
module load htslib/1.9
module load samtools

echo Aligning sample $in
echo to reference $index

bowtie2 -p 8 -x $index -U $in | samtools sort -n -O bam -o $out - 
