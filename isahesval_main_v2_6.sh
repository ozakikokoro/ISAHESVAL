#!/bin/bash

# Requirement for tools:
#Ubuntu 22.04.3 LTS
#GNU parallel 20230922
#flair 1.6.2
#samtools 1.16.1
#htslib(and bgzip, tabix) 1.16
#bcftools 1.16
#bedtools 2.31.0
#pyfaidx 0.7.2.2
#whatshap 1.7

#jvarkit f576b50de
# Please ensure that java -jar ${jvarkit_jar}, for example, invokes jvarkit program
jvarkit_jar=/your_dir/jvarkit.jar
#openjdk 17.0.8.1 2023-08-24
#OpenJDK Runtime Environment (build 17.0.8.1+1-Ubuntu-0ubuntu122.04)
#OpenJDK 64-Bit Server VM (build 17.0.8.1+1-Ubuntu-0ubuntu122.04, mixed mode, sharing)

# picard-slim 2.27.4
# pysam (I added comment every line where pysam is used, for debugging)

minimum_py=/your_dir/minimum.py
read_length_py=/your_dir/test_reconsider2.py

# This setting is for calculating read length
# Keep this setting as it is
read_length_step1_file=read_length_step1.tsv

lockfile_dir=/your_dir/lockfiles/
mkdir -p ${lockfile_dir}

# dir for nanopore .bam file
sample_class=allele
input_dir=/your_dir/
input_file=nanopore.sort.bam
# please place .bam.bai in the same folder

cd  ${input_dir}

# set directory name
dir_name="out"

# get current date in YYYYMMDD format
date_string=$(date '+%Y%m%d')

# directory pattern
pattern="${dir_name}_${date_string}"

# check if any directories exist with the current date
existing_dirs=( $pattern* )

if [ -e "${existing_dirs[0]}" ]; then
  # if directories exist, find the highest numbered suffix
  highest_suffix=0
  for dir in "${existing_dirs[@]}"; do
    suffix=${dir##*_} # get the part after the last underscore
    if [[ $suffix =~ ^[0-9]+$ ]] && [ "$suffix" -gt "$highest_suffix" ]; then # if it's a number
      highest_suffix=$suffix
    fi
  done

  # increment the highest suffix found
  new_suffix=$((highest_suffix+1))
  new_dir_name="${pattern}_${new_suffix}_v2_6"
else
  # if no directories exist, start with "_1" for consistency
  new_dir_name="${pattern}_1_v2_6"
fi

# create the new directory and check if mkdir command succeeded
if mkdir "$new_dir_name"; then
  echo "Directory $new_dir_name created successfully."
else
  echo "Error: Failed to create the directory."
  exit 1
fi

reference=/your_dir/GRCh38_no_alt_analysis_set.fasta
# Caution! please make a .dict file by the command like picard CreateSequenceDictionary -R GRCh38_no_alt_analysis_set.fasta
# Also, you have to map long-reads using the same reference for the input file (see instruction in github)

faidx ${reference} -i chromsizes > sizes.genome
# calculate the genome sizes from the reference(fasta)
# , later used for calculating gene-by-gene slopped intervals
# to be applied for removing bilateral outside-mapped chimerric reads.

splicingvar_table=/your_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.split.spliceai.filter05.converted.csv
# A table of splicing variants and gene symbols, which are selected by spliceAI evaluation
# Comma delimited
# Example: chr1,12345678,A,C,ATAD3A,12345600,12345700
# start,end are the coordinates used for samtools view, to restrict reads onto
# splicing variant and its affecting bases (see details in compare.py).

# compare.py, one of the 3 associated python scripts during preparation of the input table, assesses in intermediate4.txt by
# spliceAI score (acceptor increase, decrease, donor increase, decrease)
# its coordinate changers (like (+)1, 0, -43, -3).
# compare the four coordinates (only for the site with 0.5 or higher) and the original variant coordinate
# to make a region (1-based) start and end by taking the minimum and maximal coordinate.
# Finally this script add the two columns start(14th column) end (15th column) to the original imput file(13 columns)

# Usually prepared from spliceAI called .vcf file. and then filtered with 3 python scripts.
# This time, the file is created as follows:
# 1. DRAGEN-variant called and snpEff annotated vcf was used. Variants were further selected for heterozygous variant and PASS by DRAGEN filter
# 2. bcftools norm splitted
# 3. grep "splice" selected
# 4. SpliceAI called (annotated vcf as output)
# 5. converted into the above format (.sh plus 3 .py scripts)

snv_wgs=/your_dir/your_sample.hard-filtered_snpEff.ann.PASS_het.vcf.gz
# vcf "PASS" by DRAGEN. Also limited to "heterozygous" GT.
# This snv_wgs is used for searching candidate allele-informative SNVs

gene_symbol_coordinate=/your_dir/martquery_allgenes_4.sort.bed
# A table of chr coordinate of all the genes
# chr start end gene_symbol, tab-delimited file

geneid_genesymbol_table=/your_dir/gencode.v37_geneid_genesymbol.tsv
# geneid"/t"gene_symbol table for converting gene_symbol to geneid (for flair input)

spvar_list=${splicingvar_table}
# splicing variant list, which is applied above

threshold=10
# threshold number of long-reads for allele-informative SNVs, which is applicable to assign condition branch

genome=${reference}
annotation_gtf=/your_dir/gencode.v37.annotation.gtf
annotation_bed=/your_dir/gene_model/gencode.v37.annotation.bed
palette=/your_dir/palette01.txt
# setting for FLAIR program

# From this line on, the program begins....
# (No need to modify for each run, just modify for fine tuning)
#
#

# spvar funcion definition
spvar_sep ()
{

local list=(${1//,/ })
local chr=`echo ${list[0]}`
local coordinate=`echo ${list[1]}`
local ref=`echo ${list[2]}`
local alt=`echo ${list[3]}`
local gene_symbol=`echo ${list[4]}`
local start=`echo ${list[5]}`
local end=`echo ${list[6]}`
local spvar_id=${chr}_${coordinate}_${ref}_${alt}_${gene_symbol}

echo "----------"
echo "chr is $chr"
echo "coordinate is $coordinate"
echo "ref is $ref"
echo "alt is $alt"
echo "gene_symbol is ${gene_symbol}"
echo "start (left most predicted aberrant splicing site) is $start"
echo "end (right most predicted aberrant splicing site) is $end"
echo "..."
# To change vcf(short read WGS) file into .csv file, to check for alleles in the gene provided by splicing var.
# This subsection limit the potential allele informative variants to the region defined by gene_symbol.bed
awk -F'\t' -v gene="${gene_symbol}" -v OFS='\t' '$4 == gene' ${2} > ${spvar_id}_gene_symbol.bed
bcftools view -R ${spvar_id}_gene_symbol.bed -o ${spvar_id}_snv_wgs_pre.vcf.gz -O z ${3}
tabix ${spvar_id}_snv_wgs_pre.vcf.gz
wait

# To retrieve geneid (gencode)
local geneid=$(awk -F'\t' -v gene="${gene_symbol}" -v OFS='\t' '$2 == gene {print $1}' ${4})

# To take chr, start, end of the gene, for later using in samtools view step
local line_gene=$(head -n 1 ${spvar_id}_gene_symbol.bed)
local column1=$(echo ${line_gene} | cut -f1)
local column2=$(echo ${line_gene} | cut -f2)
local column3=$(echo ${line_gene} | cut -f3)

# To define intervals on both outsides of the gene, to filter out mapping reads
bedtools slop -i ${spvar_id}_gene_symbol.bed -g sizes.genome -l 100 -r 0 > ${spvar_id}_L1.bed
bedtools slop -i ${spvar_id}_gene_symbol.bed -g sizes.genome -l 100000 -r 0 > ${spvar_id}_L2.bed
bedtools slop -i ${spvar_id}_gene_symbol.bed -g sizes.genome -r 100 -l 0 > ${spvar_id}_R1.bed
bedtools slop -i ${spvar_id}_gene_symbol.bed -g sizes.genome -r 100000 -l 0 > ${spvar_id}_R2.bed
wait
bedtools subtract -a ${spvar_id}_L2.bed -b ${spvar_id}_L1.bed > ${spvar_id}_L3.bed
bedtools subtract -a ${spvar_id}_R2.bed -b ${spvar_id}_R1.bed > ${spvar_id}_R3.bed
wait
rm ${spvar_id}_L1.bed ${spvar_id}_L2.bed ${spvar_id}_R1.bed ${spvar_id}_R2.bed
local column1_l=$(cut -f1 ${spvar_id}_L3.bed)
local column2_l=$(cut -f2 ${spvar_id}_L3.bed)
local column3_l=$(cut -f3 ${spvar_id}_L3.bed)
local column1_r=$(cut -f1 ${spvar_id}_R3.bed)
local column2_r=$(cut -f2 ${spvar_id}_R3.bed)
local column3_r=$(cut -f3 ${spvar_id}_R3.bed)
local outside_l_region=${column1_l}:${column2_l}-${column3_l}
local outside_r_region=${column1_r}:${column2_r}-${column3_r}
wait
rm ${spvar_id}_L3.bed ${spvar_id}_R3.bed

# subsection to select hetero 0/1 0|1 1|0 genotype variants in short read WGS
bcftools query \
-f '%CHROM,%POS,%REF,%ALT\n' \
-i '(GT="0/1" || GT="1/0" || GT="0|1" || GT="1|0") && strlen(REF)==1 && strlen(ALT)==1' \
${spvar_id}_snv_wgs_pre.vcf.gz > ${spvar_id}_snv_wgs_temp1.csv

# subsection to select hetero 1/2 1|2 2|1 genotype variants in short read WGS
bcftools query \
-f '%CHROM,%POS,%ALT{0},%ALT{1}\n' \
-i 'strlen(REF)=1 && strlen(ALT[0])=1 && strlen(ALT[1])=1 && (N_ALT=2)' \
${spvar_id}_snv_wgs_pre.vcf.gz > ${spvar_id}_snv_wgs_temp2.csv

if [ ! -s "${spvar_id}_snv_wgs_temp1.csv" ] && [ ! -s "${spvar_id}_snv_wgs_temp2.csv" ]; then
        echo "Error: ${spvar_id}_snv_wgs_temp1.csv and ${spvar_id}_snv_wgs_temp2.csv did not exist or are both empty. Skipping this splicing variant..."
        rm ${spvar_id}_snv_wgs_pre.vcf.gz ${spvar_id}_snv_wgs_pre.vcf.gz.tbi ${spvar_id}_snv_wgs_temp1.csv ${spvar_id}_snv_wgs_temp2.csv \
        ${spvar_id}_gene_symbol.bed
        exit 0
        echo "Proof-of-exit0-broken"
fi
wait
cat ${spvar_id}_snv_wgs_temp1.csv ${spvar_id}_snv_wgs_temp2.csv > ${spvar_id}_snv_wgs_temp3.csv
rm ${spvar_id}_snv_wgs_pre.vcf.gz ${spvar_id}_snv_wgs_pre.vcf.gz.tbi ${spvar_id}_snv_wgs_temp1.csv ${spvar_id}_snv_wgs_temp2.csv

local SNP_list=${spvar_id}_snv_wgs_temp3.csv

# statistics info, removing empty lines and count the described lines (if 0, then output letter "0" anyway)
#grep -v '^$' ${SNP_list} | awk 'END {print NR}' >> ${7}${9}/alleleSNV_nums.txt
# : without locking
{
flock 3
grep -v '^$' ${SNP_list} | awk 'END {print NR}' >> ${7}${9}/alleleSNV_nums.txt
} 3>${12}lockfile_3

# A step to greatly reduce calculation time
# To focus on long-reads spannning on each splicing variant and cognate region
samtools view -@ 4 -b -h ${5}_filtered_pre ${chr}:${start}-${end} > ${5}_filtered_pre2_${spvar_id}
wait
sleep 5s
samtools index ${5}_filtered_pre2_${spvar_id}

# The following 3 steps were added in v.1.8 script, to exclude reads mapping to bilateral 100bp to 100kb
# outside of the gene boudaries. This treatment is to remove chimeric reads that also maps to secondary gene.
# step 1 to select alignments which maps to outside l and r 100kb regions
samtools view -@ 4 -b -h ${5}_filtered_pre2_${spvar_id} ${outside_r_region} > ${spvar_id}_rightregion.bam
samtools sort -@ 4 -m 2G ${spvar_id}_rightregion.bam > ${spvar_id}_rightregion.sort.bam
samtools index ${spvar_id}_rightregion.sort.bam

samtools view -@ 4 -b -h ${5}_filtered_pre2_${spvar_id} ${outside_l_region} > ${spvar_id}_leftregion.bam
samtools sort -@ 4 -m 2G ${spvar_id}_leftregion.bam > ${spvar_id}_leftregion.sort.bam
samtools index ${spvar_id}_leftregion.sort.bam

# step2 sam2tsv to collect outside-mapper readnames(QNAME) into text file.
java -jar ${11} sam2tsv \
-R ${6} \
${7}${spvar_id}_leftregion.sort.bam | grep -v "^#"| awk -F "\t" 'BEGIN {OFS = FS = "\t" }{ if( $6 ~ /[ACGT]/ ){ print $1}}' | uniq > ${spvar_id}_sam2tsv_left_onlyReadname.tsv

java -jar ${11} sam2tsv \
-R ${6} \
${7}${spvar_id}_rightregion.sort.bam | grep -v "^#"| awk -F "\t" 'BEGIN {OFS = FS = "\t" }{ if( $6 ~ /[ACGT]/ ){ print $1}}' | uniq > ${spvar_id}_sam2tsv_right_onlyReadname.tsv
wait
cat ${spvar_id}_sam2tsv_left_onlyReadname.tsv ${spvar_id}_sam2tsv_right_onlyReadname.tsv | uniq > ${spvar_id}_sam2tsv_bilateral_onlyReadname.tsv
rm ${spvar_id}_sam2tsv_right_onlyReadname.tsv
rm ${spvar_id}_sam2tsv_left_onlyReadname.tsv
rm ${spvar_id}_leftregion.sort.bam*
rm ${spvar_id}_rightregion.sort.bam*
rm ${spvar_id}_leftregion.bam
rm ${spvar_id}_rightregion.bam

# step3 picard to exclude the readnames which map left/right outside regions of the gene.
if [ -s "${spvar_id}_sam2tsv_bilateral_onlyReadname.tsv" ]; then
# If "${spvar_id}_sam2tsv_bilateral_onlyReadname.tsv" contains 1 or more readnames...

picard FilterSamReads \
I=${7}${5}_filtered_pre2_${spvar_id} \
O=bifilter.not_sort.${5}_${spvar_id} \
READ_LIST_FILE=${spvar_id}_sam2tsv_bilateral_onlyReadname.tsv \
FILTER=excludeReadList

samtools sort -@ 4 -m 2G bifilter.not_sort.${5}_${spvar_id} > ${5}_filtered_${spvar_id}
samtools index ${5}_filtered_${spvar_id}
rm bifilter.not_sort.${5}_${spvar_id}

else
samtools sort -@ 4 -m 2G ${7}${5}_filtered_pre2_${spvar_id} > ${5}_filtered_${spvar_id}
samtools index ${5}_filtered_${spvar_id}
fi

rm ${spvar_id}_sam2tsv_bilateral_onlyReadname.tsv
rm ${5}_filtered_pre2_${spvar_id}
rm ${5}_filtered_pre2_${spvar_id}.bai

# statistics info
#samtools view -c ${5}_filtered_${spvar_id} >> ${7}${9}/filtered_read_nums.txt
# : without lock
{
flock 4
samtools view -c ${5}_filtered_${spvar_id} >> ${7}${9}/filtered_read_nums.txt
} 4>${12}lockfile_4

# 2nd order loop for each candidate allele-informative SNV (snv_wgs_temp3.csv)
echo "----------"
echo "Evaluating candidate allele-informative SNVs:"
echo "Number of lines in ${SNP_list} is..."
awk 'END {print NR}' ${SNP_list}
echo "..."

# Process for each candidate allele-informative SNV
while read j; do
local list_2=(${j//,/ })
local chr_2=`echo ${list_2[0]}`
local coordinate_2=`echo ${list_2[1]}`
local ref_2=`echo ${list_2[2]}`
local alt_2=`echo ${list_2[3]}`
local alleleSNV=${chr_2}_${coordinate_2}_${ref_2}_${alt_2}

echo "----------"
echo "A candidate allele-informative SNV:"
echo "chr_2 is $chr_2"
echo "coordinate_2 is $coordinate_2"
echo "ref_2 (Allele 1) is $ref_2"
echo "alt_2 (Allele 2) is $alt_2"

# To convert sam to tsv for each candidate allele-informative SNV site (chr, coordinate)
# Retrieving reads which 1) originally span the splicing site and its cognate regions
# and also 2) covers the candidate allele-informative SNV site.
java -jar ${11} sam2tsv \
-R ${6} \
${7}${5}_filtered_${spvar_id} | grep -v "^#"| grep "$chr_2" | \
awk -F "\t" -v foo="${chr_2}" -v bar="${coordinate_2}" 'BEGIN {OFS = FS = "\t" }{ if($4 == foo && $8 == bar && $6 ~ /[ACGT]/ && $10 == "M" ){ print $1}}' \
> ${spvar_id}_${alleleSNV}_sam2tsv_onlyReadname.tsv
wait

if [ -s "${spvar_id}_${alleleSNV}_sam2tsv_onlyReadname.tsv" ]; then
# If sam2tsv_onlyReadname.tsv contains 1 or more readnames...
# From here, select reads from the file .._filtered
# based on readnames with indels removed above

picard FilterSamReads \
I=${7}${5}_filtered_${spvar_id} \
O=selected.not_sort.${5}_${spvar_id}_${alleleSNV} \
READ_LIST_FILE=${spvar_id}_${alleleSNV}_sam2tsv_onlyReadname.tsv \
FILTER=includeReadList
wait
samtools sort -@ 4 -m 2G selected.not_sort.${5}_${spvar_id}_${alleleSNV} > selected.${5}_${spvar_id}_${alleleSNV}
wait

samtools index selected.${5}_${spvar_id}_${alleleSNV}
rm selected.not_sort.${5}_${spvar_id}_${alleleSNV}
rm ${spvar_id}_${alleleSNV}_sam2tsv_onlyReadname.tsv

# You have .bam with indels removed
# Next, treat with samjdk again after BQ>=10 filter (change base to N with low qual)

java -jar ${11} samjdk \
-e "final byte[] quals = record.getBaseQualities();final byte[] bases = record.getReadBases();if(quals==SAMRecord.NULL_QUALS || bases==SAMRecord.NULL_SEQUENCE) return record;for(int i=0;i< quals.length;i++) {if (htsjdk.samtools.util.SequenceUtil.isNoCall(bases[i]) || (int)quals[i]  >= 10) continue;bases[i]=(byte)'N';}record.setReadBases(bases);return record;" \
selected.${5}_${spvar_id}_${alleleSNV} > ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5%.bam}.sam
wait
samtools view -@ 4 -b -h -o ${spvar_id}_${alleleSNV}_BQ10N_filtered.not_sort.${5} ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5%.bam}.sam
wait
samtools sort -@ 4 -m 2G ${spvar_id}_${alleleSNV}_BQ10N_filtered.not_sort.${5} > ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5}
wait
samtools index ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5}
wait

rm ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5%.bam}.sam
rm ${spvar_id}_${alleleSNV}_BQ10N_filtered.not_sort.${5}
rm selected.${5}_${spvar_id}_${alleleSNV} selected.${5}_${spvar_id}_${alleleSNV}.bai

# Will here consider reads without N at allele-informative SNV site
# and retrieve ref, alt QNAME
java -jar ${11} sam2tsv \
-R ${6} \
${7}${spvar_id}_${alleleSNV}_BQ10N_filtered.${5} | grep -v "^#"| grep "$chr_2" | \
awk -F "\t" -v foo="${chr_2}" -v bar="${coordinate_2}" -v refallele="$ref_2" \
'BEGIN {OFS = FS = "\t" }{ if($4 == foo && $8 == bar && $6 == refallele && $10 == "M" ){ print $1}}' \
> ${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.tsv

java -jar ${11} sam2tsv \
-R ${6} \
${7}${spvar_id}_${alleleSNV}_BQ10N_filtered.${5} | grep -v "^#"| grep "$chr_2" | \
awk -F "\t" -v foo="${chr_2}" -v bar="${coordinate_2}" -v altallele="$alt_2" \
'BEGIN {OFS = FS = "\t" }{ if($4 == foo && $8 == bar && $6 == altallele && $10 == "M" ){ print $1}}' \
> ${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.tsv

#v2_5_1_var: moving the following line in the last of the loop for allele SNV
#rm ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5} ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5}.bai

# Preparing for condition branch: the candiate allele-informative SNV has "threshold" value (e.g. 2 or 10)  or more reads for ref and alt alleles.
# If so, retrieving reads(bam) from nanopore .bam_filtered
# by using ref_2, alt_2 only QNAME(readname)

sort ${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.tsv | uniq > ${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.sort_uniq.tsv
local ref_num=$(grep -v '^$' ${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.sort_uniq.tsv | awk 'END {print NR}')
rm ${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.sort_uniq.tsv

sort ${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.tsv | uniq > ${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.sort_uniq.tsv
local alt_num=$(grep -v '^$' ${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.sort_uniq.tsv | awk 'END {print NR}')
rm ${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.sort_uniq.tsv

local total_num=$((${ref_num}+${alt_num}))

# v1.9.2 update
# If the ref and alt bases are totally supported by summed 10 or more reads(readNames),
# then the candidate allele informative SNV is OK(＝successful allele informative SNV)：
# (Ex. chr1:23456789　- splicing SNV covering long reads 10 reads (ref 2, alt 8 at allele-inf.SNV site))

# If the candidate allele informative SNV is OK, then we write it in an intermediate file
# (temp_OK_SNV.tsv).

# Also we write readNames in 2nd intermediate file(temp_OK_readname.tsv), which can be used
# to use for a conditional branch and also to make .bam file.

# To check if the haplotype block is covering the splicing variant and its cognate region,
# the 3rd intermediate file is created:temp_OK_readname_coordinate.tsv

#if [ ${ref_num} -ge ${8} ] && [ ${alt_num} -ge ${8} ]; then
if [ "${total_num}" -ge "${8}" ]; then
echo -e "${chr_2}\t${coordinate_2}" >> ${spvar_id}_temp_OK_SNV.tsv
cat ${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.tsv \
    ${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.tsv \
 >> ${spvar_id}_temp_OK_readname.tsv

cat ${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.tsv | while read line_1
do
echo -e "${line_1}\t${coordinate_2}" >> ${spvar_id}_temp_OK_readname_coordinate.tsv
done

cat ${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.tsv | while read line_2
do
echo -e "${line_2}\t${coordinate_2}" >> ${spvar_id}_temp_OK_readname_coordinate.tsv
done

# statistics info
echo "${total_num}" >> ${7}${9}/${spvar_id}_filtered2_per_alleleSNV_read_nums.txt

# ref_2 allele
picard FilterSamReads \
I=${7}${5}_filtered_${spvar_id} \
O=selected.not_sort.${5}_${spvar_id} \
READ_LIST_FILE=${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.tsv \
FILTER=includeReadList
wait

samtools sort -@ 4 selected.not_sort.${5}_${spvar_id} \
> ./${9}/${spvar_id}_selected.ref_2.${chr_2}_${coordinate_2}_${ref_2}.${5}
wait
samtools index ./${9}/${spvar_id}_selected.ref_2.${chr_2}_${coordinate_2}_${ref_2}.${5}
wait
rm selected.not_sort.${5}_${spvar_id}

# alt_2 allele
picard FilterSamReads \
I=${7}${5}_filtered_${spvar_id} \
O=selected.not_sort.${5}_${spvar_id} \
READ_LIST_FILE=${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.tsv \
FILTER=includeReadList
wait

samtools sort -@ 4 selected.not_sort.${5}_${spvar_id} \
> ./${9}/${spvar_id}_selected.alt_2.${chr_2}_${coordinate_2}_${alt_2}.${5}
wait
samtools index ./${9}/${spvar_id}_selected.alt_2.${chr_2}_${coordinate_2}_${alt_2}.${5}
wait
rm selected.not_sort.${5}_${spvar_id}

    echo "..."
    echo "This candidate SNV, ${chr_2}, ${coordinate_2}, allele 1 ${ref_2}, allele 2 ${alt_2}"
    echo "passed, and will be considered further, because total number of ref+alt ${ref_num},${alt_num} reads ${total_num} were >=${8}."
    echo "Saved allele-splitted reads which originally also cover the splicing site in 2 .bam files."

    # v2.0.0
    # Getting phaseset info of allele determining SNV (which separation of long-read is based on) if possible (in the rare cases where splicing variant itself
    # can determine the alleles of long-read of cDNA).
    if [ "${chr}" = "${chr_2}" ] && [ "${coordinate}" = "${coordinate_2}" ] && [ "${alt}" = "${alt_2}" ]; then
    local ps_allele_det_SNV=spvar_is_on_2nd_allele
    elif [ "${chr}" = "${chr_2}" ] && [ "${coordinate}" = "${coordinate_2}" ] && [ "${alt}" = "${ref_2}" ]; then
    local ps_allele_det_SNV=spvar_is_on_1st_allele
    else
    local ps_allele_det_SNV=NA
    fi
    # Writing the unique BamPair (pair of the bam files) to the BamPair_filepath table (for FLAIR).
    echo "Writing fileplaces to the BamPair_filepath.tsv"
    #echo "${spvar_id}	./${9}/${spvar_id}_selected.ref_2.${chr_2}_${coordinate_2}_${ref_2}.${5}	./${9}/${spvar_id}_selected.alt_2.${chr_2}_${coordinate_2}_${alt_2}.${5}	${geneid}${chr_2}	${coordinate_2} ${ref_2}	${alt_2}	${gene_symbol}	${chr}	${coordinate}	${ref}	${alt}	${ref_num}	${alt_num}	${ps_allele_det_SNV}" >> BamPair_filepath.tsv
    # : without locking
    {
    flock 5
    echo "${spvar_id}	./${9}/${spvar_id}_selected.ref_2.${chr_2}_${coordinate_2}_${ref_2}.${5}	./${9}/${spvar_id}_selected.alt_2.${chr_2}_${coordinate_2}_${alt_2}.${5}	${geneid}	${chr_2}	${coordinate_2}	${ref_2}	${alt_2}	${gene_symbol}	${chr}	${coordinate}	${ref}	${alt}	${ref_num}	${alt_num}	${ps_allele_det_SNV}" >> BamPair_filepath.tsv
    } 5>${12}lockfile_5

    #v2_5_1_var feature
    {
    flock 9
    # This python uses pysam
    python ${13} ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5} read_length_step1.tsv ${spvar_id}_${alleleSNV} step1_pass
    } 9>${12}lockfile_9
else
    echo "..."
    echo "This candidate SNV, ${chr_2}, ${coordinate_2}, allele 1 ${ref_2}, allele 2 ${alt_2} for ${spvar_id} did not pass"
    echo "and will not be considered further, because number of ref+alt ${ref_num},${alt_num} reads ${total_num} was <${8}."

    # statistics info
    echo "This candidate SNV, ${chr_2}, ${coordinate_2}, allele 1 ${ref_2}, allele 2 ${alt_2} for ${spvar_id} had only 1 to n (<${8}) reads"\
     >> ${7}${9}/${spvar_id}_failed_alleleSNV_insufficient_reads.txt

    #v2_5_1_var
    {
    flock 9
    # This python uses pysam
    python ${13} ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5} read_length_step1.tsv ${spvar_id}_${alleleSNV} step1_fail
    } 9>${12}lockfile_9
fi

rm ${spvar_id}_${alleleSNV}_sam2tsv_ref_2_${ref_2}_onlyReadname.tsv \
   ${spvar_id}_${alleleSNV}_sam2tsv_alt_2_${alt_2}_onlyReadname.tsv

#v2_5_1_var feature
rm ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5} ${spvar_id}_${alleleSNV}_BQ10N_filtered.${5}.bai

else
    echo "..."
    echo "This candidate SNV, ${chr_2}, ${coordinate_2}, allele 1 ${ref_2}, allele 2 ${alt_2}"
    echo "had no covering read, and was skipped."
    rm ${spvar_id}_${alleleSNV}_sam2tsv_onlyReadname.tsv

    # statistics info
    echo "This candidate SNV, ${chr_2}, ${coordinate_2}, allele 1 ${ref_2}, allele 2 ${alt_2} for ${spvar_id} had 0 read"\
    >> ${7}${9}/${spvar_id}_failed_alleleSNV_zero_read.txt
fi

done < ${SNP_list}
# This [done] is the end of per-candidate allele informative SNV

echo "----------"
echo "Splicing variant: ${chr} ${coordinate} ${ref} ${alt} for ${gene_symbol} is assessed for covering long-reads,"
echo " and also for allele-informative SNPs covered by the long-reads."
echo "..."
echo "Total number of allele informative SNVs with ${8} or more long-reads on both alleles is.."
    if [ -s "${spvar_id}_temp_OK_SNV.tsv" ]; then
    local allele_informative_SNV_num=$(awk 'END {print NR}' ${spvar_id}_temp_OK_SNV.tsv)
    echo "${allele_informative_SNV_num}"
    echo "..."
    #echo "${allele_informative_SNV_num}" >> ${7}${9}/allele_informative_SNVs_num.txt
    # : without locking
    {
    flock 6
    echo "${allele_informative_SNV_num}" >> ${7}${9}/allele_informative_SNVs_num.txt
    } 6>${12}lockfile_6

    else
    echo "..it's 0, that is, no allele informative SNVs were found for this variant."
    echo "..."

    # statistics info
    #echo "0" >> ${7}${9}/allele_informative_SNVs_num.txt
    # : without locking
    {
    flock 6
    echo "0" >> ${7}${9}/allele_informative_SNVs_num.txt
    } 6>${12}lockfile_6

    fi

    # statistics info
    local sum_f2=$(cat ${7}${9}/${spvar_id}_filtered2_per_alleleSNV_read_nums.txt | awk '{ sum_f2 += $1 } END { print sum_f2 }')
    local count_f2=$(cat ${7}${9}/${spvar_id}_filtered2_per_alleleSNV_read_nums.txt | wc -l)


    #echo ${sum_f2} >> ${7}${9}/sum_read_counts_filtered2_alleleSNV.txt
    #echo ${count_f2} >> ${7}${9}/counts_alleleSNV.txt
    #grep -v '^$' ${7}${9}/${spvar_id}_failed_alleleSNV_insufficient_reads.txt | awk 'END {print NR}' >> ${7}${9}/failed_alleleSNV_insufficient_reads_num.txt
    #grep -v '^$' ${7}${9}/${spvar_id}_failed_alleleSNV_zero_read.txt | awk 'END {print NR}' >> ${7}${9}/failed_alleleSNV_zero_read_num.txt
    # : without locking
    {
    flock 7
    echo ${sum_f2} >> ${7}${9}/sum_read_counts_filtered2_alleleSNV.txt
    echo ${count_f2} >> ${7}${9}/counts_alleleSNV.txt
    # statistics info, removing empty lines and count the line number
    grep -v '^$' ${7}${9}/${spvar_id}_failed_alleleSNV_insufficient_reads.txt | awk 'END {print NR}' >> ${7}${9}/failed_alleleSNV_insufficient_reads_num.txt
    grep -v '^$' ${7}${9}/${spvar_id}_failed_alleleSNV_zero_read.txt | awk 'END {print NR}' >> ${7}${9}/failed_alleleSNV_zero_read_num.txt
    } 7>${12}lockfile_7
    wait
    rm ${7}${9}/${spvar_id}_filtered2_per_alleleSNV_read_nums.txt \
    ${7}${9}/${spvar_id}_filtered2_per_alleleSNV_read_nums.txt \
    ${7}${9}/${spvar_id}_failed_alleleSNV_insufficient_reads.txt \
    ${7}${9}/${spvar_id}_failed_alleleSNV_zero_read.txt

if [ -s "${spvar_id}_temp_OK_SNV.tsv" ] && \
[ -s "${spvar_id}_temp_OK_readname_coordinate.tsv" ] && \
[ -s "${spvar_id}_temp_OK_readname.tsv" ]; then
    # v2.2.2 modification, to run the latter part safely when the allele informative SNVs were found
    # and remove unnecessary files smartly.
    # v2.3.0 or later does not treat conditional branch, which was needed when selecting allele SNVs which
    # encompasses spvar coordinate

    # For a conditional branch, we prepare variables
    #sort -k1 ${spvar_id}_temp_OK_readname_coordinate.tsv \
    #> ${spvar_id}_sorted_file
    #wait
    #uniq -f1 -d ${spvar_id}_sorted_file \
    #> ${spvar_id}_dup_lines

    # Here, requirement for whatshap application is, 1) presence of at least two allele-informative SNVs
    # covered by one long-read,
    # "2) the minimum and maximal coordinate of the allele-informative SNVs encompass: = disabled"
    # "the splicing variant coordinate: = disabled"

    sort ${spvar_id}_temp_OK_readname.tsv \
    > ${spvar_id}_temp_OK_readname2.tsv
    wait
    if [ "$(uniq -d -c ${spvar_id}_temp_OK_readname2.tsv | wc -l)" -gt 0 ]; then
        #  if we have at least one duplicate in the temp_OK_readname.tsv

        # :<<'#COMMENT_OUT'
        #while read line; do
        #  readname=$(echo $line | awk '{print $1}')
        #  min=$(cat ${spvar_id}_sorted_file | awk -v readname=${readname} '$1==readname{print $2}' | sort -n | head -n 1)
        #  max=$(cat ${spvar_id}_sorted_file | awk -v readname=${readname} '$1==readname{print $2}' | sort -n | tail -n 1)

          # check if the coordinate is within the min and max values, at least once in all the reads for this splicing variant.
        #  if [ "${min}" -le "${coordinate}" ] && [ "${coordinate}" -le "${max}" ]; then
        #      echo "Coordinate ${coordinate} is within the range [${min}, ${max}] for the read ${readname}"
        #      condition_branch01=yes
        #  else
        #    :
        #  fi
        #done < ${spvar_id}_dup_lines
        # COMMENT_OUT

        #if [ "${condition_branch01}" == yes ]; then
            echo "..."
            echo "For this splicing variant site, whatshap will be used to make a haplotype block"
            echo "and assign long-reads into each haplotype."
            echo "..."
            bcftools view -T <(cut -f 1,2 ${spvar_id}_temp_OK_SNV.tsv) \
            -O z ${3} > ${spvar_id}_alleleinformative_variants_for_whatshap.vcf.gz
            wait
            tabix ${spvar_id}_alleleinformative_variants_for_whatshap.vcf.gz

            # From here, whatshap is used for haplotyping and splitting all the long-reads.
            # Used data are: $reference, gene_symbol.bed(splicing variant-defined gene region),
            # and ${5}_filtered (long-reads covering splicing variant and its cognate region).
            local col1 col2 col3 col4
            IFS=$'\t' read -r col1 col2 col3 col4 < ${spvar_id}_gene_symbol.bed
            local region="${col1}:${col2}-${col3}"
            echo "$region"
            local input_nanopore_reads=${5}_filtered_${spvar_id}
            local input_vcf=${spvar_id}_alleleinformative_variants_for_whatshap.vcf.gz

            whatshap phase \
            -o ${spvar_id}_phased.vcf.gz \
            --reference=${6} \
            --ignore-read-groups \
            --tag=PS \
            ${input_vcf} \
            ${input_nanopore_reads}
            wait
            tabix ${spvar_id}_phased.vcf.gz
            wait

            whatshap haplotag \
            -o ${spvar_id}_haplotagged.bam \
            --reference ${6} \
            --ignore-read-groups \
            --output-haplotag-list ${spvar_id}_haplotag.list.gz \
            --regions ${region} \
            ${spvar_id}_phased.vcf.gz \
            ${input_nanopore_reads}
            wait
            samtools index ${spvar_id}_haplotagged.bam
            wait

            whatshap split \
            --output-h1 ./${9}/${spvar_id}_h1.bam \
            --output-h2 ./${9}/${spvar_id}_h2.bam \
            --only-largest-block \
            ${spvar_id}_haplotagged.bam \
            ${spvar_id}_haplotag.list.gz

            samtools index ./${9}/${spvar_id}_h1.bam
            samtools index ./${9}/${spvar_id}_h2.bam

            #v2_5_1_var feature
            samtools merge ./${9}/${spvar_id}_merged_h1h2.bam ./${9}/${spvar_id}_h1.bam ./${9}/${spvar_id}_h2.bam
            samtools sort -@ 4 -m 2G ./${9}/${spvar_id}_merged_h1h2.bam > ./${9}/${spvar_id}_merged_h1h2_sort.bam
            samtools index ./${9}/${spvar_id}_merged_h1h2_sort.bam
            {
            flock 9
            # Here, python uses pysam
            python ${13} ./${9}/${spvar_id}_merged_h1h2_sort.bam read_length_step1.tsv ${spvar_id}_._._h1_h2 step1_haplotyped
            } 9>${12}lockfile_9

            rm ${spvar_id}_alleleinformative_variants_for_whatshap.vcf.gz*
            mv ${spvar_id}_phased.vcf.gz ./${9}/${spvar_id}_phased.vcf.gz
            mv ${spvar_id}_phased.vcf.gz.tbi ./${9}/${spvar_id}_phased.vcf.gz.tbi
            # bgzip -c -d ./${9}/${spvar_id}_phased.vcf.gz > ./${9}/${spvar_id}_phased.vcf
            # PStrial script changed here
            #rm ${chr}_${coordinate}_${ref}_${alt}_phased.vcf.gz*

            rm ${spvar_id}_haplotag.list.gz \
               ${spvar_id}_haplotagged.bam*
            rm ./${9}/${spvar_id}_merged_h1h2.bam \
               ./${9}/${spvar_id}_merged_h1h2_sort.bam*

            # v2.0.0 feature
            bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t[ %GT]\n' \
            -r ${chr}:${coordinate} \
            -o ./${9}/${spvar_id}_phased.query-dissociated.tsv \
            ./${9}/${spvar_id}_phased.vcf.gz

                if [ -s "./${9}/${spvar_id}_phased.query-dissociated.tsv" ]; then
                  # this command creates the tsv file of a single line like:
                  # chr12   34567890       G       A        1|0
                    while read line_phase; do
                        local phasing_spvar_allele_1=$(echo ${line_phase} | awk '{print $3}')
                        local phasing_spvar_allele_2=$(echo ${line_phase} | awk '{print $4}')
                        local phasing_spvar_GT=$(echo ${line_phase} | awk '{print $5}')
                        if [ "${alt}" = "${phasing_spvar_allele_1}" ] && [ "${phasing_spvar_GT}" = "0|1" ]; then
                            local ps_allele_det_SNV=spvar_is_on_1st_haplotype
                        elif [ "${alt}" = "${phasing_spvar_allele_1}" ] && [ "${phasing_spvar_GT}" = "1|0" ]; then
                            local ps_allele_det_SNV=spvar_is_on_2nd_haplotype
                        elif [ "${alt}" = "${phasing_spvar_allele_2}" ] && [ "${phasing_spvar_GT}" = "0|1" ]; then
                            local ps_allele_det_SNV=spvar_is_on_2nd_haplotype
                        elif [ "${alt}" = "${phasing_spvar_allele_2}" ] && [ "${phasing_spvar_GT}" = "1|0" ]; then
                            local ps_allele_det_SNV=spvar_is_on_1st_haplotype
                        else
                            local ps_allele_det_SNV=error_see_variant_or_haplotype
                        fi
                    done < ./${9}/${spvar_id}_phased.query-dissociated.tsv
                else   # The case in which spvar is not on exons or, rarely, indel that is removed from candidate allele SNVs in this pipeline
                   local ps_allele_det_SNV=haplotype_unavailable_to_spvar_nonexonic_or_indel
                fi
                wait
                rm ./${9}/${spvar_id}_phased.query-dissociated.tsv \
                ./${9}/${spvar_id}_phased.vcf.gz.tbi
            # Writing the BamPair (pair of the bam files) to the BamPair_filepath table (for FLAIR).
    {
    flock 8
            echo "${spvar_id}	./${9}/${spvar_id}_h1.bam	./${9}/${spvar_id}_h2.bam	\
${geneid}	.	.	h1	h2	${gene_symbol}	\
${chr}	${coordinate}	${ref}	${alt}	${ref_num}	${alt_num}	${ps_allele_det_SNV}" \
            >> BamPair_filepath.tsv
    } 8>${12}lockfile_8

        # This block was needed for conditional branch if-close.
        #else
        # In this script, nothing will be reported for splicing variants which have 2 or more allele-informative SNVs
        # and have at least one haplotype block, which failed to encompass the splicing variant within the region.
        #echo "..."
        #echo "This splicing variant had 2 or more allele-informative SNVs, but did not have a haplotype-block \
#containing the splicing variant within."
        #echo "..."
        #fi

    else
    # This "else" is for the case in which "splicing variant had just one allele-informative SNV".
    # This script does not here describe script for splicing variants without 2 allele-informative SNVs
    # because every allele-informative SNVs are described above already.
    echo "..."
    echo "This splicing variant did not have 2 or more allele-informative SNVs."
    echo "..."
    fi

    # removing bam file especially mapping to region of interest
    # instructed in the starting part of the 1st loop
    rm ${5}_filtered_${spvar_id} \
    ${5}_filtered_${spvar_id}.bai

    # removing the splicing site specific interemediate files for each cycle
    rm ${spvar_id}_temp_OK_SNV.tsv \
       ${spvar_id}_temp_OK_readname.tsv \
       ${spvar_id}_temp_OK_readname2.tsv \
       ${spvar_id}_temp_OK_readname_coordinate.tsv \
       ${spvar_id}_gene_symbol.bed
       #${spvar_id}_sorted_file
       #${spvar_id}_dup_lines

    # removing the candidate allele-specific SNVs specifically listed for the splicing site
    rm ./${spvar_id}_snv_wgs_temp3.csv

else
# This case of "else" means "No surviving SNV in these candidate allele informative SNVs".
# Let us remove bam file especially mapping to region of interest
# instructed in the starting part of the 1st loop
rm ${5}_filtered_${spvar_id} \
   ${5}_filtered_${spvar_id}.bai
# removing the splicing site specific interemediate files for each cycle
rm ${spvar_id}_gene_symbol.bed
rm ./${spvar_id}_snv_wgs_temp3.csv
fi

}

flair_para ()
{
echo "working on ${1}..."
local PN=$(echo ${1} | awk '{print $1}')
local geneid=$(echo ${1} | awk '{print $4}')
local chr_2=$(echo ${1} | awk '{print $5}')
local coordinate_2=$(echo ${1} | awk '{print $6}')
local ref_2=$(echo ${1} | awk '{print $7}')
local alt_2=$(echo ${1} | awk '{print $8}')
local gene_symbol=$(echo ${1} | awk '{print $9}')
local chr=$(echo ${1} | awk '{print $10}')
local coordinate=$(echo ${1} | awk '{print $11}')
local ref=$(echo ${1} | awk '{print $12}')
local alt=$(echo ${1} | awk '{print $13}')
local ref_num=$(echo ${1} | awk '{print $14}')
local alt_num=$(echo ${1} | awk '{print $15}')

mkdir ${2}/${PN}_${chr_2}_${coordinate_2}_${ref_2}_${alt_2}
local output_dir=${2}/${PN}_${chr_2}_${coordinate_2}_${ref_2}_${alt_2}
local input_bam_bamtofastq1=$(echo ${1} | awk '{print $2}')
bedtools bamtofastq -i ${input_bam_bamtofastq1} -fq ${output_dir}/allele1.fastq
local input_bam_bamtofastq2=$(echo ${1} | awk '{print $3}')
bedtools bamtofastq -i ${input_bam_bamtofastq2} -fq ${output_dir}/allele2.fastq

# version2.0.0 feature
local ps=$(echo ${1} | awk '{print $16}')

# v2_5_1_var feature
local read_length=$(echo ${1} | awk '{print $17}')

echo "${3}	${ref_2}	1	./${output_dir}/allele1.fastq" > ./${output_dir}/manifest.tsv
echo "${3}	${alt_2}	2	./${output_dir}/allele2.fastq" >> ./${output_dir}/manifest.tsv

for m in `ls ./${output_dir}/allele?.fastq`
do

local sampledata_path=${m%.fastq}

flair align  \
--genome ${4} \
--reads ${m} \
--output ${sampledata_path}_flair.aligned \
--junction_bed ${5} \
-t 14
wait
sleep 0.1s

flair correct \
-q ${sampledata_path}_flair.aligned.bed \
-g ${4} \
-f ${6} \
--nvrna \
-t 14 \
--output ${sampledata_path}_flair \
--print_check
wait

done
# resulting file is ${sample_name}_flair_all_corrected.bed

cat ./${output_dir}/*_flair_all_corrected.bed | sort > ./${output_dir}/concat_flair_all_corrected.sort.bed
cat ./${output_dir}/allele?.fastq > ./${output_dir}/concat.fastq


flair collapse \
-q ./${output_dir}/concat_flair_all_corrected.sort.bed \
-g ${4} \
-f ${6} \
--reads ./${output_dir}/concat.fastq \
-t 14 \
--support 4 \
--no_gtf_end_adjustment \
--output ./${output_dir}/flair
wait
sleep 0.1s

flair quantify \
--reads_manifest ./${output_dir}/manifest.tsv \
--isoforms ./${output_dir}/flair.isoforms.fa \
-t 14 \
--output ./${output_dir}/flair

# this quantify step generates [flair.tpm.tsv] file and counts.matrix file[flair].
mv ./${output_dir}/flair ./${output_dir}/flair.counts_matrix.tsv
# This command renames simple counts.matrix file, which will be used further.
# (flair.tpm.tsv will not be used further, because only reads mapping to this region
# were used.)

plot_isoform_usage \
./${output_dir}/flair.isoforms.bed \
./${output_dir}/flair.counts_matrix.tsv \
${geneid} \
-o ./${output_dir}/flair_${geneid} \
--palette ${7}

diff_iso_usage \
./${output_dir}/flair.counts_matrix.tsv \
allele_${ref_2}_1 \
allele_${alt_2}_2 \
./${output_dir}/diff_iso_usage_result.txt
wait

rm ./${output_dir}/allele?.fastq \
   ./${output_dir}/concat.fastq

# if the p-value(3rd column) of the diff_iso_usage_result.txt is "NA", skip the line
while read preline; do
 if [ $(echo $preline | awk '{print $3}') == "NA" ]; then
    :
 else
 echo $preline >> ./${output_dir}/diff_iso_usage_result_filter.txt
 fi
done < ./${output_dir}/diff_iso_usage_result.txt

if [ -s "${output_dir}/diff_iso_usage_result_filter.txt" ]; then
    local minimum_p=$(python ${10} ./${output_dir}/diff_iso_usage_result_filter.txt)
   # local corrected_minimum_p_spvar=$(echo "${13}*${minimum_p}" | bc)
    local corrected_minimum_p_spvar=$(awk -v a="${13}" -v b="${minimum_p}" 'BEGIN{ printf "%.10e", a * b }')
    local corrected_minimum_p_bampair=$(awk -v a="${14}" -v b="${minimum_p}" 'BEGIN{ printf "%.10e", a * b }')
   # local corrected_minimum_p_bampair=$(echo "${14}*${minimum_p}" | bc)
    if [ `echo "$(printf "%.20f" $minimum_p) > 0.05" | bc` -eq 1 ] || [ `echo "$(printf "%.20f" $minimum_p) == 0.05" | bc` -eq 1 ]; then
        rm -r ./${output_dir}
        if [ "${ref_2}" = "h1" ] && [ "${alt_2}" = "h2" ]; then
            rm ./${2}/${PN}_h1.bam* ./${2}/${PN}_h2.bam*
            rm ./${2}/${PN}_phased.vcf.gz
        else
            rm ./${2}/${PN}_selected.ref_2.${chr_2}_${coordinate_2}_${ref_2}.${8}*
            rm ./${2}/${PN}_selected.alt_2.${chr_2}_${coordinate_2}_${alt_2}.${8}*
        fi
    #v2_5_1_var
        {
        flock 9
        awk -F'\t' -v key="${PN}_${chr_2}_${coordinate_2}_${ref_2}_${alt_2}" -v OFS='\t' -v cmp0="${minimum_p}" -v cmp1="${corrected_minimum_p_spvar}" -v cmp2="${corrected_minimum_p_bampair}" '$1 == key {$4 = "step2_p>=0.05"; $5 = cmp0; $6 = cmp1; $7 = cmp2; print}' ${12} >> read_length_step2.tsv
        } 9>${11}lockfile_9
    else
    #v2_5_1_var
        {
        flock 9
        awk -F'\t' -v key="${PN}_${chr_2}_${coordinate_2}_${ref_2}_${alt_2}" -v OFS='\t' -v cmp0="${minimum_p}" -v cmp1="${corrected_minimum_p_spvar}" -v cmp2="${corrected_minimum_p_bampair}" '$1 == key {$4 = "step2_p<0.05"; $5 = cmp0; $6 = cmp1; $7 = cmp2; print}' ${12} >> read_length_step2.tsv
        } 9>${11}lockfile_9
        if [ `echo "$(printf "%.20f" $corrected_minimum_p_spvar) > 0.05" | bc` -eq 1 ] || [ `echo "$(printf "%.20f" $corrected_minimum_p_spvar) == 0.05" | bc` -eq 1 ]; then
        :
        else
        {
        flock 9
        awk -F'\t' -v key="${PN}_${chr_2}_${coordinate_2}_${ref_2}_${alt_2}" -v OFS='\t' -v cmp0="${minimum_p}" -v cmp1="${corrected_minimum_p_spvar}" -v cmp2="${corrected_minimum_p_bampair}" '$1 == key {$4 = "step2_spvar_corrected_p<0.05"; $5 = cmp0; $6 = cmp1; $7 = cmp2; print}' ${12} >> read_length_step2.tsv
        } 9>${11}lockfile_9
        {
        flock 11
        echo "${PN}	${gene_symbol}	${geneid}	${chr}	${coordinate}	\
${ref}	${alt}	${chr_2}	${coordinate_2}	${ref_2}	${alt_2}	\
${ref_num}	${alt_num}	${minimum_p}	${corrected_minimum_p_spvar}	${corrected_minimum_p_bampair}	${ps}" >> PN_minpvalue005_spvar_corrected.tsv
        } 11>${11}lockfile_11
            if [ `echo "$(printf "%.20f" $corrected_minimum_p_bampair) > 0.05" | bc` -eq 1 ] || [ `echo "$(printf "%.20f" $corrected_minimum_p_bampair) == 0.05" | bc` -eq 1 ]; then
            :
            else
            {
            flock 9
            awk -F'\t' -v key="${PN}_${chr_2}_${coordinate_2}_${ref_2}_${alt_2}" -v OFS='\t' -v cmp0="${minimum_p}" -v cmp1="${corrected_minimum_p_spvar}" -v cmp2="${corrected_minimum_p_bampair}" '$1 == key {$4 = "step2_bampair_corrected_p<0.05"; $5 = cmp0; $6 = cmp1; $7 = cmp2; print}' ${12} >> read_length_step2.tsv
            } 9>${11}lockfile_9
            {
            flock 12
            echo "${PN}	${gene_symbol}	${geneid}	${chr}	${coordinate}	\
${ref}	${alt}	${chr_2}	${coordinate_2}	${ref_2}	${alt_2}	\
${ref_num}	${alt_num}	${minimum_p}	${corrected_minimum_p_spvar}	${corrected_minimum_p_bampair}	${ps}" >> PN_minpvalue005_bampair_corrected.tsv
            } 12>${11}lockfile_12
            fi
        fi
        {
        flock 4
        echo "${PN}	${gene_symbol}	${geneid}	${chr}	${coordinate}	\
${ref}	${alt}	${chr_2}	${coordinate_2}	${ref_2}	${alt_2}	\
${ref_num}	${alt_num}	${minimum_p}	${corrected_minimum_p_spvar}	${corrected_minimum_p_bampair}	${ps}" >> PN_minpvalue005.tsv
        } 4>${11}lockfile_4
    fi
        {
        flock 5
        echo "${PN}	${gene_symbol}	${geneid}	${chr}	${coordinate}	\
${ref}	${alt}	${chr_2}	${coordinate_2}	${ref_2}	${alt_2}	\
${ref_num}	${alt_num}	${minimum_p}	${corrected_minimum_p_spvar}	${corrected_minimum_p_bampair}	${ps}" >> PN_minpvalue.tsv
        } 5>${11}lockfile_5
    else
        # In case only those isoforms with "NA" exist for the allele pair, they should be removed;
        rm -r ./${output_dir}
        if [ "${ref_2}" = "h1" ] && [ "${alt_2}" = "h2" ]; then
            rm ./${2}/${PN}_h1.bam* ./${2}/${PN}_h2.bam*
        else
            rm ./${2}/${PN}_selected.ref_2.${chr_2}_${coordinate_2}_${ref_2}.${8}*
            rm ./${2}/${PN}_selected.alt_2.${chr_2}_${coordinate_2}_${alt_2}.${8}*
        fi
fi

# removing err_tmp_....txt, which is sometimes created during flair processes.
rm -f ./err_tmp_*

}

export -f spvar_sep flair_para
#↑ required for parallel tool

cd ${input_dir}
echo "Processing ${input_file} ....."

#stat output-1
stat_file=${input_dir}${new_dir_name}/Run_report.tsv
spvar_num=$(awk 'END {print NR}' ${spvar_list})
echo "# <Run report>" > ${stat_file}
echo "---------<Basic settings>----------" >> ${stat_file}
echo "# The script file itself:	${0}" >> ${stat_file}
echo "# Input_dir:	${input_dir}" >> ${stat_file}
echo "# Input_file (nanopore long-read bam file):	${input_file}" >> ${stat_file}
echo "# splicingvar_table (splicing variants examined):	${splicingvar_table}" >> ${stat_file}
echo "# (short read) SNV file (vcf) as candidates allele informative SNVs:	${snv_wgs}" >> ${stat_file}
echo "# Referece genome:	${reference}" >> ${stat_file}
echo "# Gene_symbol coordinate .bed:	${gene_symbol_coordinate}" >> ${stat_file}
echo "# GeneID gene_symbol table:	${geneid_genesymbol_table}" >> ${stat_file}
echo "# Output_dir:	${input_dir}${new_dir_name}" >> ${stat_file}
echo "----------<flair parameters>---------" >> ${stat_file}
echo "# Annotation GTF:	${annotation_gtf}" >> ${stat_file}
echo "# Annotation bed:	${annotation_bed}" >> ${stat_file}
echo "# Palette:	${palette}" >> ${stat_file}
echo "----------<other parameters>---------" >> ${stat_file}
echo "# Threshold coverage (total number of long-reads required at allele informative SNV):	${threshold}" >> ${stat_file}
echo "--------------------" >> ${stat_file}
echo "Input number of putative splicing variants (splitted):	${spvar_num}" >> ${stat_file}

samtools view -@ 60 -F 2048 -b -h -o ${input_file}_filtered_pre ${input_file}
wait
sleep 1s
samtools index -@ 16 ${input_file}_filtered_pre
# To exclude unmapped reads entirely from nanopore reads .bam

#The part added 2023.01.18, to limit only to the reads covering each SNP coordinate

# Reading a spvar_list=splicing variant table

# Note that this spvar_list can contain SNP, small insertion, and small deletion
cp ${spvar_list} ./spvar_list_tmp.csv
mapfile -t LINES < spvar_list_tmp.csv
length_array=${#LINES[@]}
echo Number of lines in "${spvar_list}" is...
echo $length_array splicing variants:
echo "${LINES[@]}"
echo "----------"

   parallel -a ${spvar_list} -j 10 --noswap --delay 0.1 --no-run-if-empty \
   spvar_sep ::: "${gene_symbol_coordinate}" ::: "${snv_wgs}" ::: "${geneid_genesymbol_table}" ::: "${input_file}" ::: "${reference}"\
   ::: "${input_dir}" ::: "${threshold}" ::: "${new_dir_name}" ::: "${stat_file}" ::: "${jvarkit_jar}" ::: "${lockfile_dir}" ::: "${read_length_py}"

#${spvar_list} ${1}
#${gene_symbol_coordinate} ${2}
#${snv_wgs} ${3}
#${geneid_genesymbol_table} ${4}
#${input_file} ${5}
#${reference} ${6}
#${input_dir} ${7}
#${threshold} ${8}
#${new_dir_name} ${9}
#${stat_file} ${10}
#${jvarkit_jar} ${11}
#${lockfile_dir} ${12}
#${read_length_py} ${13}

if [ $? -eq 1 ]; then
    echo "One or more jobs encountered errors and were skipped for [spvar_sep]."
fi

echo "----------"
echo "Completed processing of ${input_dir}${input_file}_filtered for all the splicing variants noted."
echo "Processed"
processed_var_num=$(awk 'END {print NR}' ./spvar_list_tmp.csv)
echo "${processed_var_num} splicing variants"
echo "----------"
rm ${input_file}_filtered_pre \
   ${input_file}_filtered_pre.bai \
   ./spvar_list_tmp.csv \
   sizes.genome

# statistics info, on numbers of total candidate allele informative SNV by spvar
if [ -s "${input_dir}${new_dir_name}/alleleSNV_nums.txt" ]; then
sum_nums_calleleSNVs=$(awk '/^[0-9]+$/ { sum_nums_calleleSNVs += $1 } END { print sum_nums_calleleSNVs }' "${input_dir}${new_dir_name}/alleleSNV_nums.txt")
count_nums_calleleSNVs=$(grep -c '^[0-9]\+$' "${input_dir}${new_dir_name}/alleleSNV_nums.txt")
if [ "${count_nums_calleleSNVs}" = "0" ]; then
 echo "There are no candidate allele SNVs in this dataset. Please re-check the input." >> ${stat_file}
else
 average_nums_calleleSNVs=$(echo "scale=2; ${sum_nums_calleleSNVs} / ${count_nums_calleleSNVs}" | bc)
 echo "Sum of numbers of candidate allele informative SNVs:	${sum_nums_calleleSNVs}" >> ${stat_file}
 echo "Avg. number of candidate allele informative SNVs per splicing variant:	${average_nums_calleleSNVs}" >> ${stat_file}
fi
 rm ${input_dir}${new_dir_name}/alleleSNV_nums.txt
else
 echo "Error: There are no candidate allele SNVs in this dataset, due to lack of intermediate file [alleleSNV_nums.txt]" >> ${stat_file}
fi

# statistics info, on read counts of 1st filtered long-reads for each spvar site
if [ -s "${input_dir}${new_dir_name}/filtered_read_nums.txt" ]; then
    sum_1stfilter_longreads=$(grep -v '^$' ${input_dir}${new_dir_name}/filtered_read_nums.txt | cat | awk '/^[0-9]+$/ { sum_1stfilter_longreads += $1 } END { print sum_1stfilter_longreads }')
    count_1stfilter_longreads=$(grep -v '^$' ${input_dir}${new_dir_name}/filtered_read_nums.txt | grep -c '^[0-9]\+$')
    if [ "${count_1stfilter_longreads}" = "0" ]; then
        echo "There are no long-reads covering splicing variants in the list" >> ${stat_file}
    else
        average_1stfilter_longreads=$(echo "scale=2; ${sum_1stfilter_longreads} / ${count_1stfilter_longreads}" | bc)
        echo "Sum of splicing variant-covering long-read counts:	${sum_1stfilter_longreads}" >> ${stat_file}
        echo "Avg. counts of splicing variant-covering long-reads per splicing variant:	${average_1stfilter_longreads}" >> ${stat_file}
    fi
    rm ${input_dir}${new_dir_name}/filtered_read_nums.txt
else
    echo "Error: There are no long-reads covering splicing variants in the list, due to lack of intermediate file [filtered_read_nums.txt]" >> ${stat_file}
fi

# statistics info, on 2nd per-allele SNV filtered read counts of long-reads for each spvar site
if [ -s "${input_dir}${new_dir_name}/sum_read_counts_filtered2_alleleSNV.txt" ] && [ -s "${input_dir}${new_dir_name}/counts_alleleSNV.txt" ]; then
    sum_r=$(cat ${input_dir}${new_dir_name}/sum_read_counts_filtered2_alleleSNV.txt | awk '{ sum_r += $1 } END { print sum_r }')
    sum_c=$(cat ${input_dir}${new_dir_name}/counts_alleleSNV.txt | awk '{ sum_c += $1 } END { print sum_c }')
    if [ "${sum_c}" = "0" ]; then
        echo "There are no filtered long-reads covering allele informative SNVs." >> ${stat_file}
    else
        average_r2=$(echo "scale=2; ${sum_r} / ${sum_c}" | bc)
        echo "Sum of filtered long-read counts on all allele informative SNVs for all splicing variants:	${sum_r}" >> ${stat_file}
        echo "Avg. filtered long-read counts per allele informative SNV for all splicing variants:	${average_r2}" >> ${stat_file}
    fi
    rm ${input_dir}${new_dir_name}/sum_read_counts_filtered2_alleleSNV.txt \
    ${input_dir}${new_dir_name}/counts_alleleSNV.txt
else
    echo "Error: There are no filtered long-reads covering allele informative SNVs, due to lack of intemediate files [sum_read_counts_filtered2_alleleSNV.txt counts_alleleSNV.txt]" >> ${stat_file}
fi

# statistics info, failed (allele informative SNV counts with insufficient or zero long-reads)
if [ -s "${input_dir}${new_dir_name}/failed_alleleSNV_insufficient_reads_num.txt" ]; then
    sum_failed1=$(cat ${input_dir}${new_dir_name}/failed_alleleSNV_insufficient_reads_num.txt | awk '{ sum_failed1 += $1 } END { print sum_failed1 }')
    count_failed1=$(grep -v '^$' ${input_dir}${new_dir_name}/failed_alleleSNV_insufficient_reads_num.txt | grep -c '^[0-9]\+$')
    if [ "${count_failed1}" = "0" ]; then
        echo "There are no allele informative SNVs with insufficient long-reads (1<= n <${threshold})" >> ${stat_file}
    else
        average_failed1=$(echo "scale=2; ${sum_failed1} / ${count_failed1}" | bc)
        echo "Sum of allele informative SNV counts with insufficient long-reads:	${sum_failed1}" >> ${stat_file}
        echo "Avg. of allele informative SNV counts with insufficient long-reads per splicing variant:	${average_failed1}" >> ${stat_file}
    fi
    rm ${input_dir}${new_dir_name}/failed_alleleSNV_insufficient_reads_num.txt
else
    echo "Error: There are no allele informative SNVs with insufficient long-reads (1<= n <${threshold}) \
due to lack of intermediate file [failed_alleleSNV_insufficient_reads_num.txt]" >> ${stat_file}
fi

# statistics info
if [ -s "${input_dir}${new_dir_name}/failed_alleleSNV_zero_read_num.txt" ]; then
    sum_failed2=$(cat ${input_dir}${new_dir_name}/failed_alleleSNV_zero_read_num.txt | awk '{ sum_failed2 += $1 } END { print sum_failed2 }')
    count_failed2=$(grep -v '^$' ${input_dir}${new_dir_name}/failed_alleleSNV_zero_read_num.txt | grep -c '^[0-9]\+$')
    if [ "${count_failed2}" = "0" ]; then
        echo "There are no allele informative SNVs with no long-reads (=0)" >> ${stat_file}
    else
        average_failed2=$(echo "scale=2; ${sum_failed2} / ${count_failed2}" | bc)
        echo "Sum of allele informative SNV counts with no long-reads:	${sum_failed2}" >> ${stat_file}
        echo "Avg. of allele informative SNV counts with no long-reads per splicing variant:	${average_failed2}" >> ${stat_file}
    fi
    rm ${input_dir}${new_dir_name}/failed_alleleSNV_zero_read_num.txt
else
    echo "Error: There are no allele informative SNVs with no long-reads \
due to lack of intermediate file [failed_alleleSNV_zero_read_num.txt]" >> ${stat_file}
fi

echo "Now starting preprocessing for FLAIR..."
#line 1118

uniq_num_spvar=$(cut -f1 BamPair_filepath.tsv | sort | uniq | awk 'END {print NR}')
num_bam_pairs=$(grep -v '^$' ./BamPair_filepath.tsv | awk 'END {print NR}')
# From here, FLAIR program processes the paired .bam files
# As a part 0, two .bam files will be reversed to .fastq in advance.

    parallel -a BamPair_filepath.tsv -j 10 --noswap --delay 0.1 --no-run-if-empty \
    flair_para ::: "${new_dir_name}" ::: "${sample_class}" ::: "${genome}" ::: "${annotation_bed}"\
    ::: "${annotation_gtf}" ::: "${palette}" ::: "${input_file}" ::: "${stat_file}"\
    ::: "${minimum_py}" ::: "${lockfile_dir}" ::: "${read_length_step1_file}" ::: "${uniq_num_spvar}" ::: "${num_bam_pairs}"

#${new_dir_name} ${2}
#${sample_class} ${3}
#${genome} ${4}
#${annotation_bed} ${5}
#${annotation_gtf} ${6}
#${palette} ${7}
#${input_file} ${8}
#${stat_file} ${9}
#${minimum_py} ${10}
#${lockfile_dir} ${11}
#${read_length_step1_file} ${12}
#${uniq_num_spvar} ${13}
#${num_bam_pairs} ${14}

if [ $? -eq 1 ]; then
    echo "One or more jobs encountered errors and were skipped [flair_para]."
fi

mv BamPair_filepath.tsv ./${new_dir_name}/
mv read_length_step?.tsv ./${new_dir_name}/

# statistics info on spvar and allele SNV combinations (or spvar and whatshap called h1/h2 combinations)
haplotype_pairs=$(grep "haplotype" ./${new_dir_name}/BamPair_filepath.tsv | awk 'END {print NR}')
singlesite_alleleSNV_pairs=$(grep -v "haplotype" ./${new_dir_name}/BamPair_filepath.tsv | awk 'END {print NR}')
echo "Input to flair (BamPair_filepath.tsv statistics):" >> ${stat_file}
echo "Combinations of splicing variants and allele informative SNVs (or whatshap generated haplotypes) inspected by flair:	${num_bam_pairs}" >> ${stat_file}
echo "Unique number of splicing variants inspected by flair:	${uniq_num_spvar}" >> ${stat_file}
echo "Combinations of spicing variants and single site allele informative SNVs with ${threshold} or more long-reads coverage:	${singlesite_alleleSNV_pairs}" >> ${stat_file}
echo "Combinations of Haplotypes pairs generated based on allele informative SNVs matching the criteria:	${haplotype_pairs}" >> ${stat_file}
echo "----------"

sort PN_minpvalue005.tsv > ./${new_dir_name}/PN_minpvalue005.tsv
sort PN_minpvalue005_spvar_corrected.tsv > ./${new_dir_name}/PN_minpvalue005_spvar_corrected.tsv
sort PN_minpvalue005_bampair_corrected.tsv > ./${new_dir_name}/PN_minpvalue005_bampair_corrected.tsv
sort PN_minpvalue.tsv > ./${new_dir_name}/PN_minpvalue.tsv
wait
# statistics info on successful combination or pairs of spvar and alleleSNV(or haplotype) on 3 criteria (raw/spvar-num-corrected/bampair-num-corrected)
success_005=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005.tsv | awk 'END {print NR}')
success_005_spvar_corrected=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005_spvar_corrected.tsv | awk 'END {print NR}')
success_005_bampair_corrected=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005_bampair_corrected.tsv | awk 'END {print NR}')
success_005_single=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005.tsv | grep -v "haplotype" | awk 'END {print NR}')
success_005_single_spvar_corrected=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005_spvar_corrected.tsv | grep -v "haplotype" | awk 'END {print NR}')
success_005_single_bampair_corrected=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005_bampair_corrected.tsv | grep -v "haplotype" | awk 'END {print NR}')
success_005_ht=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005.tsv | grep "haplotype" | awk 'END {print NR}')
success_005_ht_spvar_corrected=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005_spvar_corrected.tsv | grep "haplotype" | awk 'END {print NR}')
success_005_ht_bampair_corrected=$(grep -v '^$' ./${new_dir_name}/PN_minpvalue005_bampair_corrected.tsv | grep "haplotype" | awk 'END {print NR}')

uniq_spvar_p_raw=$(cut -f1 ./${new_dir_name}/PN_minpvalue005.tsv | sort | uniq | awk 'END {print NR}')
uniq_spvar_p_spvar_corrected=$(cut -f1 ./${new_dir_name}/PN_minpvalue005_spvar_corrected.tsv | sort | uniq | awk 'END {print NR}')
uniq_spvar_p_bampair_corrected=$(cut -f1 ./${new_dir_name}/PN_minpvalue005_bampair_corrected.tsv | sort | uniq | awk 'END {print NR}')


echo "Output based on flair (PN_minpvalue005.stv statistics)" >> ${stat_file}
echo "Combinations of splicing variants associated with statistically significant (raw p<0.05) isoform changes and allele informative SNVs or haplotype:	${success_005}" >> ${stat_file}
echo "---with single site allele informative SNVs:	${success_005_single}" >> ${stat_file}
echo "---with haplotype associated:	${success_005_ht}" >> ${stat_file}
echo "---Unique number of splicing variants with single allele SNV or haplotype:	${uniq_spvar_p_raw}" >> ${stat_file}
echo "----------" >> ${stat_file}

echo "Combinations of splicing variants associated with statistically significant (spvar_corrected p<0.05) isoform changes and allele informative SNVs or haplotype:	${success_005_spvar_corrected}" >> ${stat_file}
echo "---with single site allele informative SNVs:	${success_005_single_spvar_corrected}" >> ${stat_file}
echo "---with haplotype associated:	${success_005_ht_spvar_corrected}" >> ${stat_file}
echo "---Unique number of splicing variants with single allele SNV or haplotype:	${uniq_spvar_p_spvar_corrected}" >> ${stat_file}
echo "----------" >> ${stat_file}

echo "Combinations of splicing variants associated with statistically significant (bampair_corrected p<0.05) isoform changes and allele informative SNVs or haplotype:	${success_005_bampair_corrected}" >> ${stat_file}
echo "---with single site allele informative SNVs:	${success_005_single_bampair_corrected}" >> ${stat_file}
echo "---with haplotype associated:	${success_005_ht_bampair_corrected}" >> ${stat_file}
echo "---Unique number of splicing variants with single allele SNV or haplotype:	${uniq_spvar_p_bampair_corrected}" >> ${stat_file}
echo "----------" >> ${stat_file}

rm -f PN_minpvalue005.tsv PN_minpvalue.tsv PN_minpvalue005_spvar_corrected.tsv PN_minpvalue005_bampair_corrected.tsv
rm -f ./err_tmp_*
