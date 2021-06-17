#!/usr/bin/env bash

##############################
## Readme Section
##############################

#  ~~~~~~~ OVERVIEW ~~~~~~~
#  
#  Pipeline written for the Burke Lab at Oregon State U.
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
#  The pipeline uses SGE tools configured for the OSU CGRB research computing
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
## Set the reference to call against, and a VCF of known common snps to calibrate gatk against
ref=../reference/S288C_reference_sequence_R64-1-1_20110203.fasta
bqsr_vcf=../reference/SGRP2-cerevisiae-freebayes-snps-Q30-GQ30.vcf

## Set the readgroup information file - this needs to be a CSV file with all sample names
rg_info=readgroup_metadata.csv

## Maximum number of missing samples allowed to include 
## a variant in the final output
max_nmiss=0






###############################
## Job build and submission
###############################

## Are you sure you want to run?? 
echo "(NB: This script has a Readme section.) "
echo "This will submit all the jobs. Submit? (yes/no): "
read input
if [ "$input" != "yes" ]; then
  echo "Quitting."
  exit
fi


## Map fastqs to produce BAM files, submit as array job with SGE_Array
## with -b 10 (batch size of 10 -- only run this many concurrently)
## -m 50G (start jobs on machines with 50G free RAM, kill the jobs if they try to use more than this)
## -r fastq_to_vcf_logs (put output logs here, and name the job this)
echo $samples | tr ' ' '\n' | awk -v ref="$ref" -v bqsr_vcf="$bqsr_vcf" -v rg_info="$rg_info" '{print "fastqs_to_recalbam.py --fastq_folder " $1 " --reference_fasta " ref " --bqsr_vcf " bqsr_vcf " --readgroup_info " rg_info }' \
  | SGE_Array -b 10 -m 50G -r pipeline_logs/fastq_to_recalbam_logs -P 4 -q $queue


###########################
## Run haplotypecaller in gvcf mode  on a per-sample basis
samps=$(echo $samples)
for samp in $samps; do echo gatk HaplotypeCaller -R $ref -ERC GVCF --input ${samp}/haplotype_this.bam --output ${samp}/gatk_haplotypecaller_results.g.vcf.gz; done | SGE_Array -m 16G -r pipeline_logs/call_haplotypes_logs -P 1 -q $queue --hold_names pipeline_logs/fastq_to_recalbam_logs


## Joint Calling: make GenomicsDB 
echo gatk_combine.sh $ref | SGE_Array -m 50G -r pipeline_logs/GenomicsDB_logs -q $queue --hold_names pipeline_logs/call_haplotypes_logs

## Joint Calling: GenotypeGVCFs
echo gatk GenotypeGVCFs -R $ref -V combined.g.vcf.gz -O all_samples_raw.vcf.gz | SGE_Array -m 20G -r pipeline_logs/jointcall_variants_logs -q $queue --hold_names pipeline_logs/GenomicsDB_logs


## Filter Variants: Suggested in best practices to filter SNPs & INDELs separately 
## Split the results
echo 'gatk SelectVariants --variant all_samples_raw.vcf.gz --output raw_snps.vcf.gz --select-type-to-include SNP && gatk SelectVariants --variant all_samples_raw.vcf.gz --output raw_indels.vcf.gz --select-type-to-include INDEL --select-type-to-include MIXED --exclude-non-variants' \
   | SGE_Array -m 8G -r pipeline_logs/split_variant_logs -P 4 -q $queue --hold_names pipeline_logs/jointcall_variants_logs

## Filter INDELs
echo gatk_filter_variants_indels.sh raw_indels.vcf.gz filtered_indels.vcf \
   | SGE_Array -m 8G -r pipeline_logs/filter_indels_logs -P 4 -q $queue --hold_names pipeline_logs/split_variant_logs


## Filter SNPs
echo gatk_filter_variants_snps.sh raw_snps.vcf.gz filtered_snps.vcf \
   | SGE_Array -m 8G -r pipeline_logs/filter_snps_logs -P 4 -q $queue --hold_names pipeline_logs/split_variant_logs



## SNPs
echo 'java -Xmx4g -jar ../local/bin/snpEff/snpEff.jar -s annotated_snps.html -v R64-1-1.86 filtered_snps_chrm.vcf > filtered_snps_ann.vcf' | SGE_Array -m 8G -r pipeline_logs/ann_snps_logs -P 4 -q $queue --hold_names pipeline_logs/chrom_names_logs,pipeline_logs/filter_snps_logs

## INDELs
echo 'java -Xmx4g -jar ../local/bin/snpEff/snpEff.jar -s annotated_indels.html -v R64-1-1.86 filtered_indels_chrm.vcf > filtered_indels_ann.vcf' | SGE_Array -m 8G -r pipeline_logs/ann_indels_logs -P 4 -q $queue --hold_names pipeline_logs/chrom_names_logs,pipeline_logs/filter_indels_logs


## Turn annotated variants into a table
echo 'gatk VariantsToTable -V filtered_snps_ann.vcf -F CHROM -F POS -F REF -F ALT -F HOM-VAR -F HET -F HOM-REF -F NO-CALL -F TYPE -F EVENTLENGTH -F MULTI-ALLELIC -F ANN -GF GT -GF AD -O annotated_snps.txt' | SGE_Array -m 8G -r pipeline_logs/ann_snps_table_logs -P 4 -q $queue --hold_names pipeline_logs/ann_snps_logs

echo 'gatk VariantsToTable -V filtered_indels_ann.vcf -F CHROM -F POS -F REF -F ALT -F HOM-VAR -F HET -F HOM-REF -F NO-CALL -F TYPE -F EVENTLENGTH -F MULTI-ALLELIC -F ANN -GF GT -GF AD -O annotated_indels.txt' | SGE_Array -m 8G -r pipeline_logs/ann_indels_table_logs -P 4 -q $queue --hold_names pipeline_logs/ann_indels_logs

###########################

## Turn the snps VCF into a table
echo vcf_to_table.py --vcf filtered_snps.vcf --output filtered_snps.txt --num_allow_missing $max_nmiss \
   | SGE_Array -m 4G -r pipeline_logs/vcf_to_table_snps_logs -P 4 -q $queue --hold_names pipeline_logs/filter_snps_logs



## Turn the indels VCF into a table
echo vcf_to_table.py --vcf filtered_indels.vcf --output filtered_indels.txt --num_allow_missing 0 \
   | SGE_Array -m 4G -r pipeline_logs/vcf_to_table_indels_logs -P 4 -q $queue --hold_names pipeline_logs/filter_indels_logs


## Set the permissions so others in the group can read/write them
echo "chgrp -R $group_owner . && chmod -R g+w ." \
   | SGE_Array -r pipeline_logs/fix_permissions_logs --hold_names pipeline_logs/vcf_to_table_snps_logs,pipeline_logs/vcf_to_table_indels_logs,pipeline_logs/ann_snps_table_logs,pipeline_logs/ann_indels_table_logs -q $queue -m 1G -P 1

