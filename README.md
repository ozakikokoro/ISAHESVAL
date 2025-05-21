# ISAHESVAL
## <INS>Is</INS>oform <INS>a</INS>nalysis of <INS>he</INS>terozygous putative <INS>s</INS>plicing <INS>v</INS>ariants at the <INS>a</INS>llele <INS>l</INS>evel using nanopore long-read sequencing
  
### Table of cotents  

What is ISAHESVAL ?  
Tested environment  
Input requirement  
Citation  


### What is ISAHESVAL ?  
ISAHESVAL is a bioinformatics pipeline written in bash, presumed to be implemented in linux machines. It is meant for analyzing Oxford Nanopore cDNA or direct RNAseq transcriptome (or amplicon cDNA) reads with respect to splicing changes at the allele level, in combination with variant data provied by the user. It receive input files of a variant file (.vcf) and nanopore long-read cDNA/direct RNA sequencing data (.bam) which is mapped already and separate the nanopore reads into alleles, and calculate isoform using a well-known isoform analysis tool, [FLAIR](https://github.com/BrooksLabUCSC/flair).  

### Tested environment and tools
#### Environment
A workstation with  
 AMD Ryzen Threadripper 3990X 64-Core Processor  
 Memory 256Gbytes  
 (GPU is physically installed but is not required in this pipeline)  
 
### Tools  
Ubuntu 22.04.3 LTS  
GNU parallel 20230922  
flair 1.6.2  
samtools 1.16.1  
htslib(and bgzip, tabix) 1.16  
bcftools 1.16  
bedtools 2.31.0  
pyfaidx 0.7.2.2  
whatshap 1.7  
jvarkit f576b50de  
openjdk 17.0.8.1 2023-08-24  
OpenJDK Runtime Environment (build 17.0.8.1+1-Ubuntu-0ubuntu122.04)  
OpenJDK 64-Bit Server VM (build 17.0.8.1+1-Ubuntu-0ubuntu122.04, mixed mode, sharing)  
picard-slim 2.27.4  

### Required associated files  
You have to download the following files from this site into your directory where the main script exists. The main script files refer to those associated files.  
1. reference file  
One has to create .dict file using picard based on this reference file (.fasta) like this:  
```picard CreateSequenceDictionary -R GRCh38_no_alt_analysis_set.fasta```     
3. associated python script_1
4. associated python script_2

### Input requirement  

### Citation  
_under construction_
