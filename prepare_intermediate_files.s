#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --job-name="makeintermediates"
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
#### PURPOSE: Generate intermediate files assigning each read to the  ####
####          pegRNA it maps to in addition to library statistics     ####
##########################################################################

##########################################################################
#### PREREQUISITE: RUN AFTER bowtie2_align.s. Requires output from    ####
####               both forward and reverse read.                     ####
##########################################################################

##########################################################################
#### DEPENDENCY: NONE                                                 ####
##########################################################################

##########################################################################
#### OUTPUT: Generates seven output files                             ####
####         *readmap : for each input bam, contains read name and    #### 
####                    pegRNA read maps to                           #### 
####         *tokeep : readIDs and their pegRNAs that do not have     ####
####                   recombination between pegRNA and target region #### 
####         *pegids : names of all detected pegs                     ####
####         *coverage : coverage for detected pegs                   ####
####         *recombination : observed recombination rate             ####
####         *coveragestats : number pegs >=50x, >=30x, >=1x read     ####
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
#### sbatch --export=indir=INDIR prepare_intermediate_files.s         ####
#### where                                                            ####
#### #### INDIR is full path for input directory containing bam files ####
#### #### #### The basename of INDIR will appear in file outputs.     ####
##########################################################################

module purge
module load samtools

cd $indir

samp=`basename $indir`

### Create readmap files for each bam
### Uses samtools to view only mapped reads, then pulls relevant columns
for i in *bam
do
	samtools view -F 4 $i | cut -f1,3 > ${i%%.bam}.readmap
done

### Extract mapped read ids with no recombination
### Checks that readnames and mapped pegRNA are the same
paste *readmap | awk '($1 == $3) && ($2 == $4) && ($2 != "*")' | cut -f1,2 > ${samp}_tokeep.txt

### Determine recombination rate
### Calculates total number of reads which were mapped and number mapped with no recombination
### Calculates number with recombination and reports percentage
mapped=`paste *readmap | awk '($1 == $3) && ($2 != "*") && ($4 != "*")' | wc -l`
correct=`wc -l ${samp}_tokeep.txt | cut -f1 -d' '`
recom=`echo "scale=4; 1 - ( $correct / $mapped )" | bc -l`
echo ${samp} $recom > ${samp}.recombination

### Get coverage of detected peg ids
cut -f2 ${samp}_tokeep.txt | sort | uniq -c > ${samp}_coverage.txt

### Summarize coverage
thirty=`awk '$1 >= 30' ${samp}_coverage.txt | wc -l | cut -f1 -d' '`
fifty=`awk '$1 >= 50' ${samp}_coverage.txt | wc -l | cut -f1 -d' '`
total=`wc -l ${samp}_coverage.txt | cut -f1 -d' '`
echo $samp $fifty $thirty $total > ${samp}.coveragestats

### Get list of detected pegs
awk '{print $2}' ${samp}_coverage.txt > ${samp}_pegids.txt
