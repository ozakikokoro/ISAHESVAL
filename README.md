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
 
### Tools used for testing the main script
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

### Tools used for preparation of input and associated files
snpEff 5.2-1 (or other version can be applied)  
minimap2 2_2.17 (can be another version)  
spliceAI 1.3.1


### Required associated files  
You have to download the following files from this site or other public source into your directory where the main script exists. The main script files refer to those associated files.  
1. reference genome file  
   You have to download genome data (e.g. GRCh38_no_alt_analysis_set.fasta). You have to create .dict file for the tool "picard" based on this reference file (.fasta) like this:    
```picard CreateSequenceDictionary -R GRCh38_no_alt_analysis_set.fasta```  
that should be placed in the same folder with reference file (e.g. GRCh38_no_alt_analysis_set.fasta).  
3. python script_1 (minimum.py) associated with the main script: download from [this folder](python) and place it in the same directory with the main script.  
4. python script_2 (test_reconsider2.py) associated with the main script: downlowd from [this folder](python) and place it in the same directory with the main script.  
5. Gencode files: You have to download two [Gencode](https://www.gencodegenes.org/human/release_37.html) files such as gencode.v37.annotation.bed and gencode.v37.annotation.gtf
6. gene symbol - coordinate (table): A table of chr coordinate of all the genes.  
"chr  start  end  gene_symbol" tab-delimited file. You can download from [here](genome_and_gene_model) for gencode v37.  
7. geneid(Ensemble) genesymbol table  
"geneid  gene_symbol" tab-delimited file. A table for converting gene_symbol to geneid (for flair input). You can download from [here](genome_and_gene_model) for gencode v37.  
8. python script_3 (compare.py) associated with preparation of input files in the step of splicing variant filtering ("(6) Filtering with spliceAI delta score threshold and converting .vcf into .csv file"): download from [this folder](python) and place it where you do preparation of input files.  
9. python script_4 (split.py) associated with preparation of input files in the step of splicing variant filtering ("(6) Filtering with spliceAI delta score threshold and converting .vcf into .csv file"): download from [this folder](python) and place it where you do preparation of input files.  
10. python script_5 (filter.py, filter02.py, filter08.py) associated with preparation of input files in the step of splicing variant filtering ("(6) Filtering with spliceAI delta score threshold and converting .vcf into .csv file"): download from [this folder](python) and place them where you do preparation of input files. You will use one of the three files (filter.py, filter02.py, filter08.py) for conversion/filtering depending on the wanted stringency.  

### Preparation of variant files (.vcf.gz)  
Usually, short-read whole genome sequencing data (.vcf.gz), is used. In our test, we used DRAGEN Germline (e.g. v4.0.3) for mapping and variant-calling. Only the .hard-filtered.vcf.gz and index file (.vcf.gz.tbi) are needed for the next step.  
Default setting is used, such as:  
```dragen -r reference -1 read1.fastq -2 read2.fastq --enable-variant-caller true --vc-emit-ref-confidence GVCF --vc-enable-vcf-output true --output-directory out_dir --output-file-prefix prefix --enable-map-align true --enable-sv true --enable-cnv true --output-format BAM --enable-map-align-output true --enable-bam-indexing true```  

### Preparation of input files  
1. a variant table (.csv) file with base filtering (general quality check), annotation with snpEff, basic filtering with bcftools (heterozygous variant only, split when in a compound heterozygous state, splicing region), and further annotation with spliceAI. The resulting file must be transformed into further filtering of delta score threshold (such as only splicing variants with spliceAI delta score of 0.5 or higher) by python scripts. The result of step2 ("input3=/your_file_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.vcf.gz") will be used (see below) for searching allele-informative SNVs.  
Here are example scripts for preparing such files (these script may not be necessarily optimal but I wrote to make it easy to understand and follow):  
```#!/bin/bash  
#-----(1) snpEff  
input=/your_file_dir/your_sample.hard-filtered.vcf.gz  
working_dir=$(dirname $input)  
output_pre=$(basename $input)  
output_file=${output_pre%.vcf.gz}_snpEff.ann.vcf.gz

# genome version that will be used for snpEff annotation
# which is downloaded already by a command like "snpEff download hg38"
genome_version=hg38

#conda activate snpeff (if you use snpEff in conda environment)
cd $working_dir
java -Xmx31g -jar /your_dir_for_snpEff/snpEff.jar \
-v -stats report.html $genome_version $input | bcftools view -O z -o ${working_dir}/${output_file}
#conda deactivate
  
#-----(2) bcftools basic filtering such as...  
input2=/your_file_dir/your_sample.hard-filtered_snpEff.ann.vcf.gz  
input2_short=${input2%.vcf.gz}  
bcftools view \  
-i 'FILTER=="PASS" && (GT=="0/1" || GT=="1/0" || GT=="0|1" || GT=="1|0" || GT=="1/2" || GT=="2/1" || GT=="1|2" || GT=="2|1")' \  
-O z \  
-o ${input2_short}.PASS_het.vcf.gz \  
${input2}  
  
tabix ${input2_short}.PASS_het.vcf.gz  
#-----(3) bcftools splitting (converting compound heterozygous sites into heterozygous expression)  
input3=/your_file_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.vcf.gz  
input3_short=${input3%.vcf.gz}  
bcftools norm \  
-m-any \  
-O z \  
-o ${input3_short}.split.vcf.gz \  
${input3}  
  
tabix ${input3_short}.split.vcf.gz  
#-----(4) Only variants with annotations related to "splicing"  
line=/your_file_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.split.vcf.gz  
bcftools view -h $line > header.txt  
bcftools view -H -O v $line | grep "splice" > ${line%.vcf.gz}.grep_splice.vcf_pre  
cat header.txt ${line%.vcf.gz}.grep_splice.vcf_pre > ${line%.vcf.gz}.grep_splice.vcf  
bgzip -c ${line%.vcf.gz}.grep_splice.vcf > ${line%.vcf.gz}.grep_splice.vcf.gz  
tabix ${line%.vcf.gz}.grep_splice.vcf.gz  
rm ${line%.vcf.gz}.grep_splice.vcf_pre header.txt  
  
#-----(5) SpliceAI  
input5=/your_file_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.split.grep_splice.vcf  
reference=/your_references_dir/GRCh38_no_alt_analysis_set.fasta  
  
# if using conda environment  
#conda activate spliceai131  
  
spliceai -I ${input5} \  
-O ${input5%.vcf}.spliceai.vcf \  
-R ${reference} \  
-A grch38  
  
#conda deactivate  
#-----(6) Filtering with spliceAI delta score threshold and converting .vcf into .csv file  
You need three python scripts (split.py, compare.py, and filter.py) downloaded from this site into the same folder with a shell script file. filer.py will do a standard filtering of splicing variants with spliceAI's delta score of 0.5 or higher. If you prefer to use less stringent threshold, you can replace filter.py with the filter02.py, which will keep variants with delta score of 0.2 or higher. If you rather prefer more stringent threshold, you can replace filter.py with the filter08.py, which will keep variants with delta score 0.8 or higher. Here pandas module in python must be installed for calculation. The shell script file for these three python scripts (split.py, compare.py, and filter.py), for example, contains the following scripts:  
# This is the script to convert spliceai_annotated vcf file (not bgzipped)  
# into csv file (chr, coordinate, ref, alt, gene_symbol,start,end)  
# Example: chr1,12345678,A,C,ATAD3A,12345600,12345700  
# start,end are the coordinates used for samtools view, to restrict reads onto  
# splicing variant and its affecting bases (see details in compare.py).  
  
# Because spliceAI delta score of 0.5 point (acceptor increase|descrese | donor increase | decrease) is a recommended threshold,  
# we consider 0.5 point or higher score as a standard threshold for screening.  
  
# This script requires split.py, filter.py, compare.py in the following folders.  
# When you use delta score of 0.2 as threshold, you simply replace filter.py with filter02.py  

# Setting section here  
input_file=/your_file_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.split.grep_splice.spliceai.vcf  
input_dir=${input_file%/*}/  
output_file=${input_file%.vcf}.converted.delta05.csv  
  
# From here, script begins  
splitting_python_script=/same_folder_as_this_script/split.py  
filtering_python_script=/same_folder_as_this_script/filter.py  
comparing_python_script=/same_folder_as_this_script/compare.py  
  
mkdir tmp  
cd ./tmp/  
  
grep -v "^#" ${input_file} | grep "SpliceAI" | awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, substr($8, index($8, "SpliceAI")), $9, $10}' > intermediate.txt  
awk '{if ($8 ~ /,/) {split($8, a, ","); for (i in a) {sub(/^SpliceAI=/, "", a[i]); print $1,$2,$4,$5,a[i]}} else {sub(/^SpliceAI=/, "", $8); print $1,$2,$4,$5,$8}}' intermediate.txt \  
> intermediate2.txt  
  
# Set the input and output file paths  
INPUT_FILE_split=intermediate2.txt  
OUTPUT_FILE_split=intermediate3.txt  
  
# Call the Python script with the input and output file paths as arguments  
python ${splitting_python_script} $INPUT_FILE_split $OUTPUT_FILE_split  
# This script makes a 13 columns data table.  
# chr1 159927891 TAA TA IGSF9 . . . . . . . .  
# chr1 159927891 TAA T IGSF9 0.00 0.01 0.00 0.00 -4 -7 -18 -16  
# chr1 161599629 C G FCGR2C 0.98 0.13 0.00 0.00 1 30 1 -33  
  
# Set the input and output file paths  
INPUT_FILE_filter=intermediate3.txt  
OUTPUT_FILE_filter=intermediate4.txt  
# Call the Python script with the input and output file paths as arguments  
python ${filtering_python_script} $INPUT_FILE_filter $OUTPUT_FILE_filter  
# This script makes  13 columns data table  
# filtering out spliceAI scores . . . .  
# and retrieving  at least one of the four scores >=0.5.  
#chr5 83353124 G A XRCC4 0.77 0.99 0.0 0.0 7 1 32 7  
#chr5 96900192 A G ERAP2 0.0 0.0 0.0 0.51 -37 -3 20 -3  
# Set the input and output file paths  
INPUT_FILE_compare=intermediate4.txt  
OUTPUT_FILE_compare=intermediate5.txt  
# Call the Python script with the input and output file paths as arguments  
python ${comparing_python_script} $INPUT_FILE_compare $OUTPUT_FILE_compare  
# This script makes 15 columns data table  
#chr2 240950520 A G CROCC2 0.0 0.0 0.03 0.9 -42 -4 -10 -4 240950516 240950520  
#chr3 14924000 A G FGD5 0.49 0.99 0.0 0.0 8 2 0 15 14924000 14924002  
# $1,2,3,4,5,14,15 are the core data for subsequent analysis  
  
# Use cut to extract columns 1-5, 14, and 15, using space as the delimiter  
cut -d ' ' -f 1-5,14,15 "$OUTPUT_FILE_compare" > "$output_file"  
  
# Replace space delimiter with comma = last step !  
sed -i 's/ /,/g' "$output_file"  

rm intermediate.txt intermediate2.txt intermediate3.txt intermediate4.txt intermediate5.txt  
cd ../  
rmdir tmp
```  

2. a variant file (.vcf) that is an intermediate file in the above scripts (1. base filtering -> annotation with snpEff -> basic filtering (heterozygous variant only)):  
```input3=/your_file_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.vcf.gz```  
This file will be used as a variant list to search for candidate allele-informative SNVs.
  
4. a nanopore mapped reads (.bam) that is mapped against the same version of the reference genome above. Typically, mapping with minimap2 with “-ax splice -uf -k14” option
is employed. Ensure that index file such as .bam.bai is also placed in the same folder with the .bam file. Transcript reads with any modalities (direct RNA sequencing, cDNA sequencing (amplified or not), cDNA targeted amplicon sequencing) can be an input.  
  
### Calculation time and Output folder/files  
Typical calculation time: 22.6 hrs by the default parameters (input filtered by the environment described above).  
Result files are contained within a directory which is automatically created (e.g. out_20250425_1_v2_6), as follows (e.g. excerpt, produced by analysis of GM12878 direct RNA-seq data published by [Workman et al., Nat Methods 2019](https://www.nature.com/articles/s41592-019-0617-2)):  
├── chr8_89983743_T_C_NBN_chr8_89935041_C_G (<-splicing variant, Gene Symbol, and its associated allele-informative SNV. This folder contains FLAIR result of long-reads covering the splicing variant and allele-informative SNV(s))  
│   ├── allele1_flair.aligned.bam  
│   ├── allele1_flair.aligned.bam.bai  
│   ├── allele1_flair.aligned.bed  
│   ├── allele1_flair_all_corrected.bed  
│   ├── allele1_flair_all_inconsistent.bed  
│   ├── allele2_flair.aligned.bam  
│   ├── allele2_flair.aligned.bam.bai  
│   ├── allele2_flair.aligned.bed  
│   ├── allele2_flair_all_corrected.bed  
│   ├── allele2_flair_all_inconsistent.bed  
│   ├── concat_flair_all_corrected.sort.bed  
│   ├── diff_iso_usage_result_filter.txt  
│   ├── diff_iso_usage_result.txt  
│   ├── flair.counts_matrix.tsv  
│   ├── flair_ENSG00000104320.14_isoforms.png (<-diagram of isoforms)  
│   ├── flair_ENSG00000104320.14_usage.png (<-isoform usage, colored in accordance with the above diagram of isoforms)  
│   ├── flair.isoforms.bed  
│   ├── flair.isoforms.fa  
│   ├── flair.isoforms.gtf  
│   └── manifest.tsv  
├── chr8_89983743_T_C_NBN_selected.alt_2.chr8_89935041_G.NA12878-DirectRNA.pass.dedup.mm2.sort.bam (<-allele-separated reads)  
├── chr8_89983743_T_C_NBN_selected.alt_2.chr8_89935041_G.NA12878-DirectRNA.pass.dedup.mm2.sort.bam.bai  
├── chr8_89983743_T_C_NBN_selected.ref_2.chr8_89935041_C.NA12878-DirectRNA.pass.dedup.mm2.sort.bam (<-allele-separated reads)  
├── chr8_89983743_T_C_NBN_selected.ref_2.chr8_89935041_C.NA12878-DirectRNA.pass.dedup.mm2.sort.bam.bai  
├── PN_minpvalue005_bampair_corrected.tsv (<-final table corrected by number of pairs (splicing variant and allele-informative SNV or haplotype associated))  
├── PN_minpvalue005_spvar_corrected.tsv (<-final table corrected by number of splicing variants: main result)  
├── PN_minpvalue005.tsv (<-final table showing splicing variant and associated allele-informative SNV or haplotype, but not corrected for multiple testing, filtered on the p-value 0.05 as the threshold for presence of isoform changes (FLAIR-calculated p-value)  
├── PN_minpvalue.tsv (<-final table showing splicing variant and associated allele-informative SNV or haplotype, without no filtering on the p-value of presence of isoform changes (FLAIR-calculated p-value)  
├── read_length_step1_fail_nonHLA.txt (<-experimental features for read-length calculation)  
├── read_length_step1_fail.txt  
├── read_length_step1_haplotyped_nonHLA.txt  
├── read_length_step1_haplotyped.txt  
├── read_length_step1_nonHLA.tsv  
├── read_length_step1_pass_nonHLA.txt  
├── read_length_step1_passORhaplotyped_nonHLA.txt  
├── read_length_step1_passORhaplotyped.txt  
├── read_length_step1_pass.txt  
├── read_length_step1.tsv  
├── read_length_step2_fail_nonHLA.txt  
├── read_length_step2_fail.txt  
├── read_length_step2_nonHLA.tsv  
├── read_length_step2_pass_nonHLA.txt  
├── read_length_step2_pass.txt  
├── read_length_step2.tsv  
└── Run_report.tsv (<-final summary for the calculation)  
  
For quick understanding of each run, please look at Run_report.tsv and PN_minpvalue005_spvar_corrected.tsv. Then take a look at respective folder containing FLAIR-generated details for each splicing variant - allele-informative SNV pair. Only folders related to pairs with p < 0.05 isoform changes remain in the result folder ("out_yyyymmdd_v2_6").

Let us see the example of the PN_minpvalue005_spvar_corrected.tsv:  
```
id (chr and coordinate of splicing varinat and gene symbol)	gene_symbol	Ensemble_id	chr	coordinate	ref	alt	chr_of_allele-informative_SNV	coordinate_of_allele-informative_SNV	Ref_of_allele-informative_SNV	Alt_of_allele-informative_SNV	Ref-covering_long-reads	Alt-covering_long-reads	raw_p-value	p-value_corrected_by_number_of_splicing_variants	p-value_corrected_by_number_of_splicing_variant-allele-informative_SNV_combinations	haplotype_information  
chr12_109561243_C_T_MMAB	MMAB	ENSG00000139428.12	chr12	109561243	C	T	chr12	109573424	G	T	37	34	1.48E-06	2.12E-04	1.52E-03	NA  
chr12_109561243_C_T_MMAB	MMAB	ENSG00000139428.12	chr12	109561243	C	T	chr12	109573425	C	T	42	31	7.14E-07	1.03E-04	7.38E-04	NA  
chr12_109561243_C_T_MMAB	MMAB	ENSG00000139428.12	chr12	109561243	C	T	.	.	h1	h2	42	31	7.66E-11	1.10E-08	7.91E-08	spvar_is_on_1st_haplotype  
chr19_4453238_C_T_UBXN6	UBXN6	ENSG00000167671.12	chr19	4453238	C	T	chr19	4454086	C	T	119	117	5.25E-05	7.55E-03	5.42E-02	NA  
chr21_44908184_T_C_ITGB2	ITGB2	ENSG00000160255.18	chr21	44908184	T	C	.	.	h1	h2	119	6	3.31E-09	4.77E-07	3.42E-06	haplotype_unavailable_to_spvar_nonexonic_or_indel  
chr2_162279995_C_G_IFIH1	IFIH1	ENSG00000115267.9	chr2	162279995	C	G	chr2	162267541	C	T	109	76	1.83E-38	2.64E-36	1.89E-35	NA  
chr2_162279995_C_G_IFIH1	IFIH1	ENSG00000115267.9	chr2	162279995	C	G	chr2	162272314	T	C	97	103	1.19E-49	1.72E-47	1.23E-46	NA  
chr2_162279995_C_G_IFIH1	IFIH1	ENSG00000115267.9	chr2	162279995	C	G	.	.	h1	h2	1	0	2.19E-68	3.15E-66	2.26E-65	haplotype_unavailable_to_spvar_nonexonic_or_indel
```
  
P-value is the lowest p-value for isoform changes analyzed by FLAIR.  
When h1 h2 are written in the Ref/Alt of allele-informative SNV, the this line is information on whatshap created haplotype. When haplotype is created (that is, two or more allele-informative SNVs are available for this splicing variant), long-reads covering this haplotype region are all analyzed by FLAIR.  
Haplotype information for the splicing variant are often unavailable (the last column=NA or could not be calculated because the splicing variant was intronic (non-exonic) and not on the transcript long-reads).  
When haplotype information is available, to understand the relationship among the splicing variant and allele-informative SNVs, please look into {splicing variant id}.phased.vcf.gz file. In the .phased.vcf.gz file, which is created by whatshap, GT:PS tag is recorded, where phase 0|1 or 1|0 is written for each variant, in the last column. For example, if the splicing variant and allele-informative SNV 1 have 1|0 and allele-informative SNV 2 has 0|1, then the splicing variant and allele-informative SNV1 on the same chromosome, while allele-informative SNV 2 is on the other chromosome.  
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878  
chr12   109561243       .       C       T       0       PASS    KM=7.8;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka;ANN=T|intron_variant|MODIFIER|MMAB|ENSG00000139428|transcript|ENST00000545712.6|protein_coding|6/8|c.520-139G>A|||||| GT:PS    1|0:109561243  
chr12   109561263       .       G       A       0       PASS    KM=10;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka;ANN=A|intron_variant|MODIFIER|MMAB|ENSG00000139428|transcript|ENST00000545712.6|protein_coding|6/8|c.519+157C>T||||||  GT:PS    0|1:109561243  
chr12   109573424       .       G       T       0       PASS    KM=10.2;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka;ANN=T|synonymous_variant|LOW|MMAB|ENSG00000139428|transcript|ENST00000545712.6|protein_coding|1/9|c.57C>A|p.Arg19Arg|451/4438|57/753|19/250||,T|upstream_gene_variant|MODIFIER|MVK|ENSG00000110921|transcript|ENST00000228510.7|protein_coding||c.-1399G>T|||||388|,T|intragenic_variant|MODIFIER|MVK|ENSG00000110921|gene_variant|ENSG00000110921|||n.109573424G>T||||||     GT:PS  1|0:109561243  
chr12   109573425       .       C       T       0       PASS    KM=9.94;KFP=0;KFF=0;MTD=bwa_freebayes,bwa_gatk,bwa_platypus,isaac_strelka;ANN=T|missense_variant|MODERATE|MMAB|ENSG00000139428|transcript|ENST00000545712.6|protein_coding|1/9|c.56G>A|p.Arg19His|450/4438|56/753|19/250||,T|upstream_gene_variant|MODIFIER|MVK|ENSG00000110921|transcript|ENST00000228510.7|protein_coding||c.-1398C>T|||||387|,T|intragenic_variant|MODIFIER|MVK|ENSG00000110921|gene_variant|ENSG00000110921|||n.109573425C>T||||||  GT:PS  1|0:109561243  
```
  


### Citation  
_under construction_
