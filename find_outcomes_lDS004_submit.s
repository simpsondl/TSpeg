#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --job-name="editing_outcomes"
#SBATCH --output=/Genomics/grid/users/ds65/logs/peg_outcomes/%x-%j.out

module load R/4.2.0

echo "Now processing fastq file $fq for cell line $cl, sample day $day, replicate $repl"
timestart=`date`

echo "Starting at $timestart"

Rscript ~/scripts/pegRNA/editing_pipeline.R \
--file=$fq \
--pegID=$id \
--cellline=$cl \
--day=$day \
--replicate=$repl \
--barcode=$bc \
--targetregion=$tr \
--output=$out

timeend=`date`

echo "Completed processing at $timeend"
