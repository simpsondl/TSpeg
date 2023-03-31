#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=50G
#SBATCH --job-name="makepegfq"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out

##########################################################################
#### Before running this, edit the SBATCH directives above to point   ####
#### to your own log directory and to have appropriate time and       #### 
#### memory constraints based on the amount of sequencing performed   #### 
#### and number of sequences that need to be aligned. Additionally,   ####
#### adjust the path of the script make_peg_fastas_mini.s to point to ####
#### a local copy of this pipeline.                                   ####
##########################################################################

##########################################################################
#### PURPOSE: Split original FASTQ by pegRNA each read maps to        ####
##########################################################################

##########################################################################
#### PREREQUISITE: RUN AFTER prepare_intermediate_files_driver.s      ####
##########################################################################

##########################################################################
#### OUTPUT: One trimmed and one untrimmed FASTQ containing target    ####
####         region readouts for each pegRNA in library               ####    
##########################################################################

##########################################################################
#### This script expects a directory and number of pegRNAs to process ####
#### in each job to be given as input. The input directory should     ####
#### contain the outputs resulting from prepare_intermediate_files.   #### 
#### The number of pegRNAs in each job, set by the sz variable, and   ####
#### the total number of pegRNAs in a library together determine the  #### 
#### number of "chunks" that get processed. Temp files for each chunk ####
#### will be made at this step and a job is created for each one. Use ####
#### appropriate setting for sz to avoid creating too many jobs and   ####
#### reaching limits. This script can be run after a partial run that ####
#### processes only some of the pegRNAs for a sample.                 ####
##########################################################################

##################################################################################################
#### Call script with                                                                         ####
#### sbatch --export=indir=INDIR,sz=SIZE,sampfq=SAMPLEFQ make_peg_fastas_driver.s             ####
#### where                                                                                    ####
#### #### INDIR is full path for input directory containing outputs from PREREQS              ####
#### #### SIZE is the number of pegRNAs to process in each job                                ####
#### #### #### LIBRARY_SIZE / SIZE = NUM JOBS SUBMITTED                                       ####
#### #### SAMPLEFQ is the full path for the sample FASTQ file containing the target region    ####
##################################################################################################

module purge
module load htslib/1.9

cd $indir

### Create output directories
mkdir -p peg_fastas/trimmed peg_fastas/untrimmed

### Line of last pegid
line=`ls peg_fastas/trimmed | wc -l`
### Number of lines remaining
tosplit=`tail -n +${line} ${indir%/}_pegids.txt | wc -l`
### Number of files to make
nchunks=$(( ($tosplit / $sz) + 1 ))

### Create file with pegids still to be processed
tail -n +${line} ${indir%/}_pegids.txt > tmp.txt

### Create split pegid files
for i in `seq $nchunks`
do 
	start=$(( (($i - 1) * $sz + 1) ))
	stop=$(( $i * $sz ))
	sed -n "${start},${stop}p;${stop}q" tmp.txt > tmp_${i}.txt 
done

for i in tmp_*
do
	sbatch --export=in1=${indir},in2=${sampfq},in3=${i} \
	~/scripts/pegRNA/make_peg_fastas_mini.s
done
