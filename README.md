# TSpeg
Analysis pipeline used in Cirincione and Simpson et al. ("A benchmarked, high-efficiency prime editing platform for multiplexed dropout screening") to determine editing efficiencies from included target sites in the self-targeting library. 

### Processing Pipeline
Run the scripts in the following order to process self-targeting library data. Detailed instructions are located within each file.

1. bowtie2_align.s
2. prepare_intermediate_files.s
3. make_peg_fastas_driver.s
4. find_outcomes_lDS004_driver.s

### Running Pipeline on Demo Data
This pipeline is intended to be run on a SLURM-based HPC environment to enable high-throughput processing of multiple samples and thousands of epegRNAs simultaneously. In the demo data, a selection of 3 epegRNAs from the full dataset of 2,000 will be processed.

1. Place all scripts from this repository (files ending in .s or .R) into a local directory on a SLURM HPC system.
2. Download demo fastq files (demo_r1.fq and demo_r2.fq) as well as the demo_peg_info.txt file into the same directory as all scripts. Note: Different directory structures can be accomodated provided that scripts are updated to internally use the full paths for all dependent scripts; comments in each file indicate file paths that would need to be adjusted.
3. Make a new directory to contain demo output files.
   ```
   mkdir demo_output
   ```
4. Run bowtie2_align.s on each fastq file independently. Note: The output path includes the directory made in the previous step, adjust as necessary.
   ```
   sbatch --export=in=demo_r1.fq,index=reference_files/demo_r1_reference,out=demo_output/demo_r1.bam bowtie2_align.s
   sbatch --export=in=demo_r2.fq,index=reference_files/demo_r2_reference,out=demo_output/demo_r2.bam bowtie2_align.s
   ```
5. Run prepare_intermediate_files.s, supplying the path to the output folder. 
   ```
   sbatch --export=indir=demo_output prepare_intermediate_files.s
   ```
6. Run make_peg_fastas_driver.s. This script will spawn a new task for each epegRNA (3 in the demo data). Note: The script moves into the output directory, so the sampfq argument here uses a relative path to point to the correct location; this can be adjusted to use a full path in a local copy. Additionally, while the trimb and trime options are provided for flexibility, further steps do not use the resulting trimmed fastq files as inputs.
   ```
   sbatch --export=indir=demo_output/,sz=1,sampfq=../demo_r2.fq,trimb=6,trime=6 make_peg_fastas_driver.s
   ``` 
