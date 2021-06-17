#!/usr/bin/env bash

## Dropping this in a script is mostly to avoid all of the crazy quoting etc. required by the filter strings...

if [ $# -ne 2 ]; then
  echo "This script parses an input vcf, which has been selected for SNPs, with gatk VariantFiltration to produce an output with some specific filters."
  echo "Usage: gatk_filter_variants.sh <input.vcf.gz> <output.vcf.gz>"
  exit
fi

gatk VariantFiltration --variant $1 --output $2 --filter-name LowVQCBD --filter-expression "QD < 5.0" --filter-name LowQualSNPs --filter-expression "MQ < 40.0"  --filter-name FisherStrand --filter-expression "FS > 60.0" --filter-name ReadPosRank --filter-expression "ReadPosRankSum < -8.0" --filter-name MQRankSum --filter-expression "MQRankSum < -12.5"
