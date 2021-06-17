#!/usr/bin/env python2.7
import argparse
import io
import re
import sys
import os
import os.path
import subprocess
import re

#########################################################
## Parse input, and figure out where we'll be working
parser = argparse.ArgumentParser(description = "Turns a VCF-file into an R-importable file. Keeps only bi-allelic snps marked as PASS.")
parser.add_argument('--vcf', required = True, dest = 'vcf', help = "VCF File to parse.")
parser.add_argument('--output', required = True, dest = 'output', help = "File to write (as tab-separated table).")
parser.add_argument('--num_allow_missing', required = False, dest = 'num_allow_missing', default=0, type=int, help = "Allow missing data for up to this many samples; any more will cause the variant to be skipped in the output.")

args = parser.parse_args()


fhandle = io.open(args.vcf, "rU")
outhandle = io.open(args.output, "wb")
samples = list()

for line in fhandle:
    line_stripped = line.strip()
    
    # Sample names are encoded from field 10 to N in the ID line of the VCF
    if line_stripped.startswith("#CHROM"):
        line_list = re.split(r"\s+", line_stripped)
        for i in range(9, len(line_list)):
            samples.append(line_list[i])
            
        # write the header field
        outhandle.write("Nmiss\tCHROM\tPOS\tREF\tALT")
        for sample in samples:
            outhandle.write("\talt_" + sample + "\tN_" + sample)
        outhandle.write("\n")

    # For each variant...
    elif not line_stripped.startswith("#"):
        line_list = re.split(r"\s+", line_stripped)
        
        # if it is marked as PASS
        if line_list[6] == "PASS":

            # and if it is bi-allelic
            if not re.search(r",", line_list[4]):
                
                # parse the format str (field 9) to make sure we get the AD (allele depth) and DP (total depth) fields
                # (note: this isn't the official definition of these fields, and may depend on the program producing the VCF)
                format_list = re.split(r":", line_list[8])    # e.g. ["GT", "AD", "DP", "PL", "GQ"]
                ids_to_index = dict()                         # want to build e.g. {"GT": 0, "AD": 1, "DP": 2, "PL": 3, "GQ": 4}
                for i in range(0, len(format_list)):
                    ids_to_index[format_list[i]] = i
                
                # grab the counts
                samples_to_counts = dict()                    # samples to minor allele count, total depth; e.g. {"Sample1": [4, 12], "Sample2": [3, 18], "Sample3": [2, 14]}
                nmiss = 0                                     # number of samples missing data
                for i in range(9, len(line_list)):
                    sample_name = samples[i - 9]
                    sample_i_info = line_list[i]
                    if sample_i_info == "./.":                # missing data for this sample
                        nmiss = nmiss + 1
                        samples_to_counts[sample_name] = ["NA", "NA"]
                                        
                    else:
                        sample_i_list = re.split(r":", sample_i_info)
                        dp_count = sample_i_list[ids_to_index["DP"]]
                        ad_counts = sample_i_list[ids_to_index["AD"]]
                        minor_count = re.split(r",", ad_counts)[1]
                        samples_to_counts[sample_name] = [minor_count, dp_count]


                # write it out
                if nmiss <= args.num_allow_missing:
                    chrom = line_list[0]
                    pos = line_list[1]
                    ref = line_list[3]
                    alt = line_list[4]
                    outhandle.write(str(nmiss) + "\t" + chrom + "\t" + pos + "\t" + ref + "\t" + alt)
                    for sample in samples:
                        minor_count = samples_to_counts[sample][0]
                        dp_count = samples_to_counts[sample][1]
                        outhandle.write("\t" + minor_count + "\t" + dp_count)
                    outhandle.write("\n")
                
                
                else:
                    sys.stderr.write("Warning: skipping variant with: " + str(nmiss) + " missing variants: " + line_stripped + "\n")
                    

            else:
                sys.stderr.write("Warning: skipping non-biallelic variant: " + line_stripped + "\n")


fhandle.close()    
outhandle.close()
