#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --job-name="makepegfq"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out


#### Call script with                 
#### sbatch --export=indir=INDIR,in1=PEGIDS.TXT,in2=SAMPLE.FQ.GZ,in3=DEMUXDIR

module purge
module load htslib/1.9

/Genomics/grid/users/ds65/scripts/pegRNA/make_peg_fastas_lDS003_parallel.sh $in1 $in2 $in3 $in4
