#!/usr/bin/env bash

##############################
## Readme Section
##############################

#  ~~~~~~~ OVERVIEW ~~~~~~~
#
#  This pipeline is part 2 and should be run after run_pipeline-add-samples-01.sh is finished
#  and the resulting sample folder (including the *.g.vcf file) are moved into
#  the organized_fastqs folder that contains the rest of the samples. This script 
#  starts a joint calling and runs through the end of the pipeline producing annotated 
#  variant tables. 
#
#  Pipeline written for the Burke Lab at Oregon State U. by Shawn O'Neil@CGRB@OSU, March 2018
#  
#  This pipeline runs GATK 4.X on a set of populations/samples to produce two 
#  tables of parsed and filtered variant calls: filtered_snps.txt and 
#  filtered_indels.txt. The outputs are tab-separated files of this format:
#  
#  Nmiss  CHROM    POS      REF  ALT  alt_sample1_sm  N_sample1_sm  alt_sample1_sm  N_sample2_sm  alt_sample2_sm
#  0      chr1     1294     C    T    1               21            1               38            1
#  0      chr1     1302     A    G    8               22            16              39            10
#  0      chr1     1303     G    T    3               22            10              40            9
#  0      chr1     1337     G    A    5               20            8               41            4
#  
#  Here, alt_sample1_sm is the count of alternate alleles for sample1, and 
#  N_sample1_sm is the total count for sample1. Nmiss is the number of samples 
#  with missing data (encoded as NA). 
#  
#  
#  ~~~~~~~ CAVEATS ~~~~~~~
#  
#  The pipeline uses SGE tools configured for the CGRB research computing
#  cluster, and is unlikely to work outside this environment.
#  
#  Currently the pipeline only supports *single-end* fastq or fastq.gz inputs. 
#  Supporting paired end would require re-architecting fastqs_to_vcf.py. 
#  
#  Currently the pipline hard-codes some gatk variant filtering criteria
#  in gatk_filter_variants.sh.
#  
#  It only outputs to the table *bi-allelic* variants; a warning is 
#  produced on stderr (in pipeline_logs/vcf_to_table_indels_logs/*.e*) for variants
#  that don't meet this.
#  
#  Currently the pipeline *requires* a reference VCF file for GATK base-quality-
#  score-recalibration ('bqsr'). This requirement could be removed by modifying this 
#  script.
#  
#  It's possible to have the preliminary steps (per-sample fastq -> bam)
#  use scratch space, by setting --scratch_dir in the calls to fastqs_to_vcf.py below.
#  I've left this off for now, since if many sub-jobs land on the same
#  machine (very likely) then there may not be enough scratch to 
#  accommodate them all. 
#  
#  
#  ~~~~~~~ HOW TO USE ~~~~~~~
#  
#  I'd suggest making a "project" folder with three main subfolders:
#  local/bin (containing a copy of the pipeline's scripts and binaries, including fastqs_to_vcf.py 
#  and others), organized_fastqs/ (for this pipeline file and sample 
#  folders), and reference/ (for the reference fasta and the BQSR calibration VCF). 
#  The reference .fasta and BQSR .vcf will need to be indexed in various ways,
#  but if they are not the fastqs_to_vcf.py will complain and let you know
#  what to do to index them (in pipeline_logs/fastq_to_vcf_logs/*.e*)
#  
#  The .fastq or .fastq.gz files could be soft-links to the real files.
#  (e.g. with `ln -s /global/fastq/storage/sample1.fastq.gz organized_fastqs/sample1/`)
#  
#  Example layout:
#  
#   project_name/
#  	    local/bin/{fastqs_to_vcf.py, vcf_to_table.py, gatk, etc}
#   
#  	    organized_fastqs/      
#  	        run_pipeline.sh    (this file)
#           POP_sample1/
#  	            input1.fastq.gz
#  	            input2.fastq.gz
#  	            input3.fastq.gz
#  	        POP_sample2/
#  	            input4.fastq.gz
#  	            input5.fastq.gz
# 
#  	    reference/
#  	        ref.fasta
#  	        bqsr_vcf.vcf
#  
#  
#  Next, edit the config section below to match the paths and sample names,
#  and run this script on the CGRB login machine (shell.cgrb.oregonstate.edu):
#  
#  cd project_name/organized_fastqs
#  ./run_pipeline.sh
#  
#  Note that the samples= line should match the sample folder names and nothing else,
#  for the above example we'd want something like samples=POP_*
#
#


##############################
## Config section
##############################

## add path to binaries and scripts
export PATH=


## Pattern to use to match sample folders
samples=POP_*


## Set the reference to call against, and a VCF of known common snps to calibrate gatk against
ref=../reference/S288C_reference_sequence_R64-1-1_20110203.fasta
bqsr_vcf=../reference/SGRP2-cerevisiae-freebayes-snps-Q30-GQ30.vcf

## Set the readgroup information file
rg_info=readgroup_metadata.csv

## Maximum number of missing samples allowed to include 
## a variant in the final output
max_nmiss=0

## Iteration variable (run pipeline multiple time)
iter_var="iter_2"

mkdir ${iter_var}

## SGE Queue string
queue=


## When it's all said and done, change everything to be writable by this group 
group_owner=burke_lab





###############################
## Job build and submission
###############################

## Are you sure you want to run?? 
echo "(NB: This script has a Readme section.) "
echo "This will submit all the jobs starting at joint calling. Submit? (yes/no): "
read input
if [ "$input" != "yes" ]; then
  echo "Quitting."
  exit
fi



## Joint Calling: make GenomicsDB 
echo gatk_combine.sh $ref ${iter_var} | SGE_Array -m 50G -r pipeline_logs_${iter_var}/GenomicsDB_logs -q $queue 

## Joint Calling: GenotypeGVCFs
echo gatk GenotypeGVCFs -R $ref -V ${iter_var}/combined.g.vcf.gz -O ${iter_var}/all_samples_raw.vcf.gz | SGE_Array -m 20G -r pipeline_logs_${iter_var}/jointcall_variants_logs -q $queue --hold_names pipeline_logs_${iter_var}/GenomicsDB_logs


## Filter Variants: Suggested in best practices to filter SNPs & INDELs separately 
## Split the results
echo "gatk SelectVariants --variant ${iter_var}/all_samples_raw.vcf.gz --output ${iter_var}/raw_snps.vcf.gz --select-type-to-include SNP && gatk SelectVariants --variant ${iter_var}/all_samples_raw.vcf.gz --output ${iter_var}/raw_indels.vcf.gz --select-type-to-include INDEL --select-type-to-include MIXED --exclude-non-variants" \
   | SGE_Array -m 8G -r pipeline_logs_${iter_var}/split_variant_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/jointcall_variants_logs

## Filter INDELs
echo gatk_filter_variants_indels.sh ${iter_var}/raw_indels.vcf.gz ${iter_var}/filtered_indels.vcf \
   | SGE_Array -m 8G -r pipeline_logs_${iter_var}/filter_indels_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/split_variant_logs


## Filter SNPs
echo gatk_filter_variants_snps.sh ${iter_var}/raw_snps.vcf.gz ${iter_var}/filtered_snps.vcf \
   | SGE_Array -m 8G -r pipeline_logs_${iter_var}/filter_snps_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/split_variant_logs


## Run SNPEff
## Chromosome names are different chrm1 vs. I, need to change 
echo "cat ${iter_var}/filtered_snps.vcf | sed 's/^chr1\t/I\t/' | sed 's/^chr2\t/II\t/' | sed 's/^chr3\t/III\t/' | sed 's/^chr4\t/IV\t/' | sed 's/^chr5\t/V\t/' | sed 's/^chr6\t/VI\t/' | sed 's/^chr7\t/VII\t/' | sed 's/^chr8\t/VIII\t/' | sed 's/^chr9\t/IX\t/' | sed 's/^chr10\t/X\t/' | sed 's/^chr11\t/XI\t/' | sed 's/^chr12\t/XII\t/' | sed 's/^chr13\t/XIII\t/' | sed 's/^chr14\t/XIV\t/' | sed 's/^chr15\t/XV\t/' | sed 's/^chr16\t/XVI\t/' | sed 's/^chrmito\t/Mito\t/' > ${iter_var}/filtered_snps_chrm.vcf && cat ${iter_var}/filtered_indels.vcf | sed 's/^chr1\t/I\t/' | sed 's/^chr2\t/II\t/' | sed 's/^chr3\t/II\t/' | sed 's/^chr4\t/IV\t/' | sed 's/^chr5\t/V\t/' | sed 's/^chr6\t/VI\t/' | sed 's/^chr7\t/VII\t/' | sed 's/^chr8\t/VIII\t/' | sed 's/^chr9\t/IX\t/' | sed 's/^chr10\t/X\t/' | sed 's/^chr11\t/XI\t/' | sed 's/^chr12\t/XII\t/' | sed 's/^chr13\t/XIII\t/' | sed 's/^chr14\t/XIV\t/' | sed 's/^chr15\t/XV\t/' | sed 's/^chr16\t/XVI\t/' | sed 's/^chrmito\t/Mito\t/' > ${iter_var}/filtered_indels_chrm.vcf" | SGE_Array -m 8G -r pipeline_logs_${iter_var}/chrom_names_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/filter_snps_logs,pipeline_logs_${iter_var}/filter_indels_logs


## SNPs
echo "java -Xmx4g -jar ../local/bin/snpEff/snpEff.jar -s ${iter_var}/annotated_snps.html -v R64-1-1.86 ${iter_var}/filtered_snps_chrm.vcf > ${iter_var}/filtered_snps_ann.vcf" | SGE_Array -m 8G -r pipeline_logs_${iter_var}/ann_snps_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/chrom_names_logs,pipeline_logs_${iter_var}/filter_snps_logs

## INDELs
echo "java -Xmx4g -jar ../local/bin/snpEff/snpEff.jar -s ${iter_var}/annotated_indels.html -v R64-1-1.86 ${iter_var}/filtered_indels_chrm.vcf > ${iter_var}/filtered_indels_ann.vcf" | SGE_Array -m 8G -r pipeline_logs_${iter_var}/ann_indels_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/chrom_names_logs,pipeline_logs_${iter_var}/filter_indels_logs


## Turn annotated variants into a table
echo "gatk VariantsToTable -V ${iter_var}/filtered_snps_ann.vcf -F CHROM -F POS -F REF -F ALT -F HOM-VAR -F HET -F HOM-REF -F NO-CALL -F TYPE -F EVENTLENGTH -F MULTI-ALLELIC -F ANN -GF GT -GF AD -O ${iter_var}/annotated_snps.txt" | SGE_Array -m 8G -r pipeline_logs_${iter_var}/ann_snps_table_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/ann_snps_logs

echo "gatk VariantsToTable -V ${iter_var}/filtered_indels_ann.vcf -F CHROM -F POS -F REF -F ALT -F HOM-VAR -F HET -F HOM-REF -F NO-CALL -F TYPE -F EVENTLENGTH -F MULTI-ALLELIC -F ANN -GF GT -GF AD -O ${iter_var}/annotated_indels.txt" | SGE_Array -m 8G -r pipeline_logs_${iter_var}/ann_indels_table_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/ann_indels_logs

###########################

## Turn the snps VCF into a table
echo vcf_to_table.py --vcf ${iter_var}/filtered_snps.vcf --output ${iter_var}/filtered_snps.txt --num_allow_missing $max_nmiss \
   | SGE_Array -m 4G -r pipeline_logs_${iter_var}/vcf_to_table_snps_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/filter_snps_logs



## Turn the indels VCF into a table
echo vcf_to_table.py --vcf ${iter_var}/filtered_indels.vcf --output ${iter_var}/filtered_indels.txt --num_allow_missing 0 \
   | SGE_Array -m 4G -r pipeline_logs_${iter_var}/vcf_to_table_indels_logs -P 4 -q $queue --hold_names pipeline_logs_${iter_var}/filter_indels_logs


## Set the permissions so others in the group can read/write them
echo "chgrp -R $group_owner . && chmod -R g+w ." \
   | SGE_Array -r pipeline_logs_${iter_var}/fix_permissions_logs --hold_names pipeline_logs_${iter_var}/vcf_to_table_snps_logs,pipeline_logs_${iter_var}/vcf_to_table_indels_logs,pipeline_logs_${iter_var}/ann_snps_table_logs,pipeline_logs_${iter_var}/ann_indels_table_logs -q $queue -m 1G -P 1


