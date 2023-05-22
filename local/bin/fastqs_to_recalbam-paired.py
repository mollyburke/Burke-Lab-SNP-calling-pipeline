#!/usr/bin/env python2.7
import argparse
import io
import re
import sys
import os
import os.path
import subprocess
import re
import csv

#########################################################
## Parse input, and figure out where we'll be working
parser = argparse.ArgumentParser(description = "Given a folder of .fastq or .fastq.gz files, as well as a reference genome (bwa and gatk indexed), and optionally a VCF file for GATK-BQSR recalibration (gatk indexed), produces a .vcf file in that same folder. Basically a wrapper for bwa and gatk HaplotypeCaller. NOTE: Currently only supports single-end reads.")
parser.add_argument('--fastq_folder', required = True, dest = 'fastq_folder', help = "Folder containing .fastq or .fastq.gz files; this folder will also be the 'working directory' for various steps of the bwa+gatk pipeline, unless --scratch_dir is set. The output .vcf file will show up in here as well. NOTE: currently each .fastq or .fastq.gz file will be mapped *independently*, ie, we don't support paired-end information (yet).")
parser.add_argument('--reference_fasta', required = True, dest = 'reference_fasta', help = "Reference genome fasta file. Must end in .fasta or .fa (gatk requirement).")
parser.add_argument('--bqsr_vcf', required = False, dest = 'bqsr_vcf', help = "VCF file for gatk's base-quality-score-recalibration process. May produce more accurate variant calls. Recommended by GATK best practices. Optional.")
parser.add_argument('--scratch_dir', required = False, dest = 'scratch_dir', help = "Use this directory for scratch work. Maybe try --scratch_dir $TMPDIR to use the system-defined scratch space.")
parser.add_argument('--readgroup_info', required = True, dest = 'readgroup_info'\
        , help = "csv file with the following columns: \
        Flow cell,Lane,Population,Renamed fastq,Original fastq,Sequencing facility,\
        Internal library name,Run type,# reads \
        Note: 'Renamed fastq files' MUST match names of fastq files in organized fastq \
        files folder- including file extention.")
args = parser.parse_args()

scratch_dir = args.fastq_folder
if args.scratch_dir:
    scratch_dir = args.scratch_dir 

print("Working directory is: " + scratch_dir)
print("Input/output directory is: " + args.fastq_folder)

notrailing = re.subn(r"/+$", "", args.fastq_folder, 0)[0]
samplename = os.path.basename(notrailing)
print("Sample name is: " + notrailing)

#########################################################
## Get list of all the fastq files in the folder, check that there are any
all_files = os.listdir(args.fastq_folder)

fastq_files = list()
for file in all_files:
    if re.search(r"(_R1.fastq$)|(_R1.fastq.gz$)", file):
        fastq_files.append(file)

assert len(fastq_files) > 0, "Error: the folder given, " + args.fastq_folder + ", seems to contain no _R1.fastq or _R1.fastq.gz files?"

#########################################################
## Check that the reference exists, is properly named, and is properly indexed
assert os.path.isfile(args.reference_fasta), "Error: the reference fasta, " + args.reference_fasta + ", doesn't seem to exist?"
assert re.search(r"(\.fasta$)|(\.fa$)", args.reference_fasta), "Error: the reference fasta, " + args.reference_fasta + ", must end in .fa or .fasta. (gatk requirement, sorry)."
dict_name = re.sub(r"(\.fasta$)|(\.fa$)", ".dict", args.reference_fasta)
assert os.path.isfile(dict_name), "Error: the reference fasta, " + args.reference_fasta + ", does not have an .dict index.\n Try running `gatk CreateSequenceDictionary --REFERENCE " + args.reference_fasta + "`"
assert os.path.isfile(args.reference_fasta + ".fai"), "Error: the reference fasta, " + args.reference_fasta + ", does not have an .fai index.\n Try running `samtools faidx " + args.reference_fasta + "`"
assert os.path.isfile(args.reference_fasta + ".bwt"), "Error: the reference fasta, " + args.reference_fasta + ", does not have an .bwt index.\n Try running `bwa index " + args.reference_fasta + "`"



#########################################################
## If given a VCF for BQSR, be sure it exists and is properly indexed
if args.bqsr_vcf:
    assert os.path.isfile(args.bqsr_vcf), "Error: the BQSR VCF, " + args.bqsr_vcf + ", doesn't seem to exist?"
    assert os.path.isfile(args.bqsr_vcf + ".idx"), "Error: the BQSR VCF, " + args.bqsr_vcf + ", does not have a .idx index.\n Try running `gatk IndexFeatureFile --feature-file " + args.bqsr_vcf + "`"
    

#########################################################
# Get readgroup info. Read from readgroup_metadata.csv file
if args.readgroup_info:
    assert os.path.isfile(args.readgroup_info), "Error: the read group information file, " + args.readgroup_info + ", does not seem to exist?"

rg_dic = {}
with open(args.readgroup_info,'rb') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        flowcell = row[0]
        lane = row[1]
        samp = row[2]
        fq = row[3]
        bc = row[4]
        facility = row[5]
        lib = row[6]
        ID = flowcell + ".lane-" + lane + "." + bc
        rg = "@RG\\tID:" + ID + "\\tSM:" + samp + "\\tPL:illumina\\tLB:" + lib + \
                "\\tPU:" + flowcell + "-" + bc + "." + lane
        rg_dic[fq] = rg


#########################################################
## Generate sam files by running bwa mem, keep a list of sam files generated
## read input from main input, write output to scratch_dir (even if that is still the same location)
## from here on out, we'll work in scratch_dir


samfiles = list()
for fastq in fastq_files:
    samfile = fastq + ".sam"
    fq2 = re.sub(r'_R1', '_R2', fastq)
    cmd = "bwa mem -t 4 -M -R \"" + rg_dic[fastq] + "\" " + args.reference_fasta + " "  + args.fastq_folder + "/" + fastq + " "  + args.fastq_folder + "/" + fq2 + " -o " + scratch_dir + "/" + samfile
    print(cmd)
    subprocess.check_output(cmd, shell = True)
    samfiles.append(samfile)


#########################################################
## Merge sam files, and delete originals
cmdlist = ['gatk MergeSamFiles']
for samfile in samfiles:
    cmdlist.append('--INPUT ' + scratch_dir + "/" + samfile)

cmdlist.append('--OUTPUT ' + scratch_dir + "/merged_from_sams.bam")
cmd = " ".join(cmdlist)
subprocess.check_output(cmd, shell = True)

for samfile in samfiles:
    cmd = "rm -f " + scratch_dir + "/" + samfile
    subprocess.check_output(cmd, shell = True)


#########################################################
## Mark Duplicates
## remove temps
cmd = "gatk MarkDuplicates --INPUT " + scratch_dir + "/merged_from_sams.bam --METRICS_FILE " + scratch_dir + "/duplicate_metrics.txt --OUTPUT " + scratch_dir + "/duplicates_marked.bam"
subprocess.check_output(cmd, shell = True)
cmd = "rm -f " + scratch_dir + "/{merged_from_sams.bam}"
subprocess.check_output(cmd, shell = True)


#########################################################
## Now, if they've specified a VCF for BQSR, we need to compute the adjustments
## otherwise, we'll just move the /dups_marked.bam to the file to be actually analyzed
if args.bqsr_vcf:
    cmd = "gatk BaseRecalibrator --input " + scratch_dir + "/duplicates_marked.bam --known-sites " + args.bqsr_vcf + " --output " + scratch_dir + "/recalibration_metrics.table --reference " + args.reference_fasta
    subprocess.check_output(cmd, shell = True)
    cmd = "gatk ApplyBQSR --bqsr-recal-file " + scratch_dir + "/recalibration_metrics.table --input " + scratch_dir + "/duplicates_marked.bam --output " + scratch_dir + "/haplotype_this.bam" 
    subprocess.check_output(cmd, shell = True)
    cmd = "rm -f " + scratch_dir + "/duplicates_marked.bam"
    subprocess.check_output(cmd, shell = True)
else:
    cmd = "mv " + scratch_dir + "/duplicates_marked.bam " + scratch_dir + "/haplotype_this.bam"
    subprocess.check_output(cmd, shell = True)


#########################################################
## Index that bam file
cmd = "gatk BuildBamIndex --INPUT " + scratch_dir + "/haplotype_this.bam"
subprocess.check_output(cmd, shell = True)


#########################################################
## Finally, move the results back to the main directory, if scratch_dir is not already the main directory
if scratch_dir != args.fastq_folder:
    # if there was no bqsr, this will show an error on trying to move recalibration_metrics.table, but that shouldn't hurt anything
    cmd = "mv " + scratch_dir + "/{haplotype_this.bam,haplotype_this.bam.bai,recalibration_metrics.table,duplicate_metrics.txt} " + args.fastq_folder
    subprocess.check_output(cmd, shell = True)
