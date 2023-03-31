#!/bin/bash

# Takes a directly containing two sorted bam files as input
# Generates seven output files:
### *readmap : for each bam, file containing readname and reference id read mapped to
### *tokeep : readids not showing recombination and reference id
### *pegids : identity of sequenced pegs
### *coverage : coverage for detected pegs
### *recombination : observed [strict] recombination rate
### *coveragestats : number pegs covered >=50x, >=30x, and >=1x 

module load samtools

cd $1

samp=${1##*/}

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
