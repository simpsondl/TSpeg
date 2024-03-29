#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --job-name="outcome_driver"
#SBATCH --output=%x-%j.out

#############################################
# Purpose: Script calls editing outcome     #
# pipeline for all pegRNAs in a sample.     #
# Submits one job for each pegRNA. Call     #
# this script after running make_peg_fastas #
#############################################

################################################################################################################
# Submission call:                                                                                             #
# sbatch --export=indir=DIR,peginfo=FILE,cl=STRING,day=INT,repl=INT find_outcomes_lDS004_driver.s              #
### where:                                                                                                     #
##### indir = sample directory containing peg_fastas directory output by make_peg_fastas script                #
##### peginfo = file containing 3 columns, no header: pegID, barcode, target region from meta_info excel sheet #
##### cl = cell line; this string will appear in plot titles and output tables                                 #
##### day = day sample collected, integer; this number will appear in plot titles and output tables            #
##### repl = replicate number, integer; this number will appear in plot titles and output tables               # 
################################################################################################################

timestart=`date`

echo Starting job at $timestart

cd $indir

for i in peg_fastas/untrimmed/*fq
do
	samp=${indir##*/}
	id=${i##*/}
	id=${id%.fq}
	id=${id%_"$samp"}
	bc=`grep $id $peginfo | cut -f2`
	tr=`grep $id $peginfo | cut -f3`

	out=editing_outcomes/${samp}/${id}
	mkdir -p $out

	echo Submitting job for pegRNA $id for sample $samp
	echo Expected barcode is $bc
	echo Expected target region $tr

	sbatch --export=fq=$i,id=$id,cl=$cl,day=$day,repl=$repl,bc=$bc,tr=$tr,out=$out ../find_outcomes_lDS004_submit.s
done	
	
timeend=`date`

echo Completed job at $timeend

