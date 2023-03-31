#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --job-name="psalign"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out


#### Call script with                 
#### sbatch --export=indir=WORKDIR,index=REF,in=READ.FQ,out=OUT bwa_protospacer.s

module purge
module load htslib/1.9
module load samtools

cd $indir

echo Aligning sample $in
echo to reference $index

bowtie2 -p 8 -x $index \
-U $in | samtools sort -n - > $out
