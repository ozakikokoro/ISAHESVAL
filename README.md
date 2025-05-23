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
You have to download the following files from this site into your directory where the main script exists. The main script files refer to those associated files.  
1. reference file
   You have to create .dict file for the tool "picard" based on this reference file (.fasta) like this:    
```picard CreateSequenceDictionary -R GRCh38_no_alt_analysis_set.fasta```  
that should be placed in the same folder with reference file (e.g. GRCh38_no_alt_analysis_set.fasta).  
3. associated python script_1: download from this site and place it in the same directory with the main script.  
4. associated python script_2: downlowd from this site and place it in the same directory with the main script.  
5. Gencode files: You have to download two [Gencode](https://www.gencodegenes.org/human/release_37.html) files such as gencode.v37.annotation.bed and gencode.v37.annotation.gtf
6. 

### Preparation of input files  
1. a variant table (.csv) file with base filtering (general quality check), annotation with snpEff, basic filtering with bcftools (heterozygous variant only, split when in compount heterozygous, splicing region), and further annotation with spliceAI. The resulting file must be transformed into further filtering of delta score threshold (such as only splicing variants with spliceAI delta score of 0.5 or higher) by python scripts.  
Here are example scripts for preparing such files (but may not work in your environment):  
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
#-----(5) Filtering with spliceAI delta score threshold and converting .vcf into .csv file  
You need three python scripts (split.py, compare.py, and filter.py) downloaded from this site into the same folder with a shell script file. filer.py will do a standard filtering of splicing variants with spliceAI's delta score of 0.5 or higher. If you prefer to use less stringent threshold, you can replace filter.py with the filter02.py, which will keep variants with delta score of 0.2 or higher. If you rather prefer more stringent threshold, you can replace filter.py with the filter08.py, which will keep variants with delta score 0.8 or higher. Here pandas module in python must be installed for calculation. The shell script file for these three python scripts (split.py, compare.py, and filter.py), for example, contains the following scripts:  
```# This is the script to convert spliceai_annotated vcf file (not bgzipped)  
# into csv file (chr, coordinate, ref, alt, gene_symbol,start,end)  
# Example: chr1,12345678,A,C,ATAD3A,12345600,12345700  
# start,end are the coordinates used for samtools view, to restrict reads onto  
# splicing variant and its affecting bases (see details in compare.py).  
  
# Because spliceAI score of 0.5 point (acceptor increase|descrese | donor increase | decrease) is a recommended threshold,  
# we consider 0.5 point or higher score as a standard threshold for screening.  
  
# This script requires split.py, filter.py, compare.py in the following folders.  
# When you use delta score of 0.2 as threshold, you simply replace filter.py with filter02.py  

# Setting section here  
input_file=/your_file_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.split.grep_splice.spliceai.vcf  
input_dir=${input_file%/*}/  
output_file=${input_file%.vcf}.converted.delta05.csv  
  
# From here, script begins  
splitting_python_script=/same_folter_as_this_script/split.py  
filtering_python_script=/same_folter_as_this_script/filter.py  
comparing_python_script=/same_folter_as_this_script/compare.py  
  
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
rmdir tmp```  

3. a variant file (.vcf) that is an intermediate file in 1. base filtering -> annotation with snpEff -> basic filtering (heterozygous variant only).
4. a nanopore mapped reads (.bam) that is mapped against the same version of the reference genome. Typically, mapping with minimap2 is employed. Ensure that index file such as .bam.bai is also placed in the same folder with the .bam file.


### Citation  
_under construction_
