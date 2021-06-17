#!/usr/bin/env bash

## This script is to start the joint genotyping process of all samples. 
## Must run GenomicsDB per chromosome and places everything in a GenomicsDB  

if [ $# -ne 1 ]; then
  echo "This script combines .g.vcf.gz, in preperation for joint calling, usings sample-name-map file to list all g.vcf.gz files. The GenomicsDB can then be used to produce a raw vcf.gz file."
  echo "Usage: gatk_joint_call.sh <refgenome> "
  exit
fi


ref=$1

## Make list of g.vcfs 
find . -name "*.g.vcf.gz" > cohort.sample_map.list

# Combine g.vcf files

gatk CombineGVCFs -R $ref -V cohort.sample_map.list -O combined.g.vcf.gz



