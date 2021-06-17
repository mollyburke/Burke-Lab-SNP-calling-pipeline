# Burke-Lab-SNP-calling-pipeline
Scripts for use with Pool-SEQ data to generate tables of SNP frequencies as txt files.


This pipeline was written for the Burke Lab at Oregon State U by Shawn O'Neill, Dana Gibbon, and Molly Burke.

This pipeline runs GATK 4.X on a set of populations/samples to produce two tables of parsed and filtered variant calls: filtered_snps.txt and 
filtered_indels.txt.

For users at OSU, the run_pipeline.sh script can be customized (see the readme within) and used with SGE tools configured for OSU's CGRB research computing cluster.  As written, this script is unlikely to work outside this environment.

Users outside of the OSU environment should reference the scripts in local/bin/ for details of the filtering parameters used with GATK, as well as how the VCF files generated by GATK are converted into SNP tables in txt format.
