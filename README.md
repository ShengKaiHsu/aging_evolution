This is the repository for the manuscript "The evolution of aging transcriptome in Drosophila"
Scripts for the data analyses and visualization are stored here.  

* [keyData](keyData/) contains the important phenotypic data and RNASeq counts for all subsequent statistical analyses  
* For the raw read data, they'll be uploaded to ENA under the accession ID: XXXXXX; (for now, they're stored locally on vetgrid26:/Volumes/Temp2/shengkai/age_CGE/formal_experiment/Pool_{518..521}a_*.fq.gz)  
* [trimming_and_mapping.sh](scripts/trimming_and_mapping.sh) is the script for read trimming, demultiplexing, mapping and QC.  
* [Age_formal_counts_generation.r](scripts/Age_formal_counts_generation.r) generates the gene counts from bam files.  
* [age_CGE_formal.R](scripts/age_CGE_formal.R) is the script to run all main analyses and visualization (Note the path).
* [usefulTable](usefulTable/) contains some utility files that are called in scripts/age_CGE_formal.R

