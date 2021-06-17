#!/usr/bin/env bash

## This script is to start the joint genotyping process of all samples. 
## Must run GenomicsDB per chromosome and places everything in a GenomicsDB  

if [ $# -ne 3 ]; then
  echo "This script imports .g.vcf.gz into a GenomicsDB usinga sample-name-map file to list all g.vcf.gz files. The GenomicsDB can then be used to produce a raw vcf.gz file."
  echo "Usage: gatk_joint_call.sh <outdatabase name> <refgenome> <list of chromosomes>"
  exit
fi


ref=$2
int=$3

## Make list of g.vcfs 
gvcfs=$(find . -name "*.g.vcf.gz")
samp=$(for i in $gvcfs; do echo "$i" | cut -d '/' -f 2 ; done)
paste <(printf "%s\n" "${samp[@]}") <(printf "%s\n" "${gvcfs[@]}") > cohort.sample_map


# GenomicsDBImport

while read line; do
	gatk GenomicsDBImport -R $ref --genomicsdb-workspace-path $1 --sample-name-map cohort.sample_map --reader-threads 8 --batch-size 50 -L $line
done < $int 



