#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=50G
#SBATCH --job-name="makepegfq"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out


#### Call script with                 
#### sbatch --export=indir=INDIR,sz=SIZE make_peg_fastas_lDS004_parallel_driver.s

module purge
module load htslib/1.9

cd $indir

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
	tmp="tmp.txt"
	sed -n "${start},${stop}p;${stop}q" tmp.txt > ${tmp%%.txt}_${i}.txt 
done

for i in tmp_*
do
	sbatch --export=in1=${indir},in2="2114__${indir%/}-read-3.fastq.gz",in3="/Genomics/adamsonlab/data/HTS_data/HTS_AC013/demultiplexed",in4=${i} ~/scripts/pegRNA/make_peg_fastas_lDS003_parallel_mini.s
done
