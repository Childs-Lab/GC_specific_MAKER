#!/usr/bin/python

#This script relies on two files produced from the MAKER_GC_cutoff_determination.pl
#script: _gc_content.txt and _cutoff.txt files.
#This script will produce two GFF files for high and low GC training.

#Tiffany Liu
#November 10, 2016

import argparse
import os.path
import sys
import os
import glob

parser = argparse.ArgumentParser(description="This program produces two GFF files for high and low GC training.")
parser.add_argument("--input_file_gff", help="Path to maker gff.", required = True)
parser.add_argument("--input_file_GC_content", help="Path to _gc_content.txt file.", required = True)
parser.add_argument("--input_file_GC_cutoff", help="Path to _cutoff.txt file.", required = True)
parser.add_argument("--output_file_low", help="Path to the low GC GFF file.", required = True)
parser.add_argument("--output_file_high", help="Path to the high GC GFF file.", required = True)
parser.add_argument("--genome_fasta", help="Path to the genome fasta file.", required = True)
args = parser.parse_args()

input_file_gff = args.input_file_gff
if os.path.isfile(input_file_gff) == False:
    print '\n\nThe file ' + input_file_gff + ' does not exist.\n'
    sys.exit()
input_file_GC_content = args.input_file_GC_content
if os.path.isfile(input_file_GC_content) == False:
    print '\n\nThe file ' + input_file_GC_content + ' does not exist.\n'
    sys.exit()
input_file_GC_cutoff = args.input_file_GC_cutoff
if os.path.isfile(input_file_GC_cutoff) == False:
    print '\n\nThe file ' + input_file_GC_cutoff + ' does not exist.\n'
    sys.exit()
genome_fasta = args.genome_fasta
if os.path.isfile(genome_fasta) == False:
    print '\n\nThe file ' + genome_fasta + ' does not exist.\n'
    sys.exit()
    
output_file_high = args.output_file_high
if os.path.isfile(output_file_high):
    print '\n\nThe file ' + output_file_high + ' already exists.\n'
    sys.exit()
output_file_low = args.output_file_low
if os.path.isfile(output_file_low):
    print '\n\nThe file ' + output_file_low + ' already exists.\n'
    sys.exit()

output_fh_high = open(output_file_high, 'w')
output_fh_low = open(output_file_low, 'w')

with open(input_file_GC_cutoff) as input_fh_GC_cutoff:
    for each_line in input_fh_GC_cutoff:
        cutoff_string = each_line.strip()
        cutoff_tuple = cutoff_string.split()
        low_cutoff = int(cutoff_tuple[0]) #Storing low GC cutoff value
        high_cutoff = int(cutoff_tuple[-1]) #Sotring high GC cutoff value
print "Low GC cutoff: ", low_cutoff
print "High GC cutoff: ", high_cutoff

lowGClist = []
lowGCset = set(lowGClist)
highGClist = []
highGCset = set(highGClist)
with open(input_file_GC_content) as input_fh_GC_content:
    for each_line in input_fh_GC_content:
        if each_line[0:6] == "seq_id": #Skipping header
            pass
        else:
            content_string = each_line.strip()
            content_tuple = content_string.split()
            seq_id = content_tuple[0]
            seq_gene = seq_id.split('-mRNA')
            seq_gene_id = seq_gene[0]
            seq_mRNA = seq_id.split('_CDS')
            seq_mRNA_id = seq_mRNA[0]
            length = int(content_tuple[5])
            GC = int(content_tuple[-1])
            if (length >= 800) and (length <= 2400): #Filtering by length
                if GC <= low_cutoff:
                    lowGCset.add(seq_gene_id)
                    lowGCset.add(seq_mRNA_id)
                if GC >= high_cutoff:
                    highGCset.add(seq_gene_id)
                    highGCset.add(seq_mRNA_id)
            else:
                pass #length is <800 or >2400
print "Length of low GC set: ", len(lowGCset)
print "Length of high GC set: ", len(highGCset)


lowCount = 0
highCount = 0
with open(input_file_gff) as input_fh_gff:
    for each_line in input_fh_gff:
        if "##FASTA" in each_line[0:7]:
            break #Stop when you reach the FASTA header
        else:
            if "#" in each_line[0]: #Keeping headers
                output_fh_high.write(each_line)
                output_fh_low.write(each_line)
            else:
                gff_string = each_line.strip()
                gff_tuple = gff_string.split()
                source = gff_tuple[1] #We are only interested in "." and maker
                if source == ".":
                    output_fh_high.write(each_line)
                    output_fh_low.write(each_line)
                elif source == "maker":
                    feature = gff_tuple[2]
                    description = gff_tuple[-1]
                    desc_id = description.split(";")
                    ID = desc_id[0]
                    identity = ID[3:]
                    parent = desc_id[1]
                    mRNAid = parent[7:]
                    if (feature == "gene") or (feature == "mRNA"):
                        if identity in lowGCset:
                            lowCount = lowCount + 1
                            output_fh_low.write(each_line)
                        if identity in highGCset:
                            highCount = highCount + 1
                            output_fh_high.write(each_line)                        
                    else: #feature == exon, 3'UTR, 5'UTR, or CDS
                        if mRNAid in lowGCset:
                            output_fh_low.write(each_line)
                        if mRNAid in highGCset:
                            output_fh_high.write(each_line)
print "Total instances of low GC genes and mRNA: ", lowCount
print "Total instances of high GC genes and mRNA: ", highCount

output_fh_low.write("%s\n"%("##FASTA")) #Adding fasta header
output_fh_high.write("%s\n"%("##FASTA"))

with open(genome_fasta) as genome_fasta_fh: #Adding genome fasta
    for each_line in genome_fasta_fh:
        output_fh_low.write(each_line)
        output_fh_high.write(each_line)

output_fh_low.close()
output_fh_high.close()
            
        
    
