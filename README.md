# TSpeg
Analysis pipeline used in Cirincione and Simpson et al. ("A benchmarked, high-efficiency prime editing platform for multiplexed dropout screening") to determine editing efficiencies from included target sites in the self-targeting library. 

### Processing Pipeline
Run the scripts in the following order to process self-targeting library data.

1. bowtie2_align.s
2. prepare_intermediate_files.s
3. make_peg_fastas_driver.s
4. find_outcomes_driver.s
