#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --job-name="makepegfq"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out

##########################################################################
#### BEFORE RUNNING THIS, edit the SBATCH directives above to point   ####
#### to your own log directory and to have appropriate time and       ####
#### memory constraints based on the amount of sequencing performed   ####
#### and number of sequences that need to be aligned. Additionally,   ####
#### check the DEPENDENCY section below for scripts that are called   ####
#### by this file, then adjust file paths to point to a local copy of ####
#### this pipeline and make any necessary changes to those scripts.   ####
##########################################################################

##########################################################################
#### PURPOSE: Acts as intermediate between driver that chunks work    ####
####          and actual processing script                            ####
##########################################################################

##########################################################################
#### PREREQUISITE: DO NOT USE THIS SCRIPT DIRECTLY. This script is    ####
####               called from inside make_peg_fastas_driver.s        ####
##########################################################################

##########################################################################
#### DEPENDENCY: Calls the script make_peg_fastas_parallel.sh         ####
##########################################################################

##########################################################################
#### OUTPUT: One trimmed and one untrimmed FASTQ containing target    ####
####         region readouts for each pegRNA in library               ####
##########################################################################

#############################################################################################################
#### Call script with                                                                                    ####
#### sbatch --export=sampfq=SAMPLEFQ,chnk=CHUNK,trimb=TRIM5PRIME,trime=TRIM3PRIME make_peg_fastas_mini.s ####
#### where                                                                                               ####
#### #### SAMPLEFQ is the full path for the sample FASTQ file containing the target region               ####
#### #### CHUNK is the local path to the tmp chunk file containing pegRNA IDs to process                 ####
#### #### TRIM5PRIME is the number of bases to remove from front/beginning of the read                   ####
#### #### TRIM3PRIME is the number of bases to remove from the end of the read                           ####
#############################################################################################################

module purge
module load htslib/1.9

ndx=${chnk##tmp_}
currdir=`pwd`
samp=`basename $currdir`

while read p; 
do 
	if [[ $p == '*' ]]; 
	then 
		continue; 
	else 
		grep $p ${samp}_tokeep.txt | cut -f1 > keeptmp_${ndx}; 
		~/tools/seqtk/seqtk subseq \
		$sampfq \
		keeptmp_${ndx} > peg_fastas/untrimmed/${p}_${samp}.fq;
		~/tools/seqtk/seqtk trimfq -b $trimb -e $trime peg_fastas/untrimmed/${p}_${samp}.fq \
		> peg_fastas/trimmed/${p}_${samp}_trimmed.fq
	fi; 
done < $chnk
