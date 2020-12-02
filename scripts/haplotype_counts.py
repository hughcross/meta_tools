#!/usr/bin/env python

# coding: utf-8
# # Extract Haplotype sequences from bam files and create table
import sys
import os
import pysam
from allele_functions.cigar_funcs import *
import argparse
# ## Get variables needed (use Argparse)

parser = argparse.ArgumentParser(description='count haplotypes from list of bam files')

parser.add_argument('-b', '--bam_list', dest='input',
type=str,
help='file with list of bam files')

parser.add_argument('-p', '--positions', dest='posifile',
type=str,
help='file with haplotype SNP positions on alignment')

parser.add_argument('-H', '--haplotypes', dest='haplos',
type=str,
help='tab-delimited file with haplotype name (1st col) and hap string (2nd col)')

parser.add_argument('-r', '--ref_name', dest='ref',
type=str,
help='Reference Name')

parser.add_argument('-d', '--directory', dest='directory',
type=str,
help='location of input bams [optional: default is current dir]')

parser.add_argument('-o', '--output_table', dest='out',
type=str,
help='output table name (default: haplotype_count_table.tsv)')

### argparse parameters
args = parser.parse_args()
bamFileList = args.input
hapFileName = args.haplos
posFileName = args.posifile
refName = args.ref
outfile = args.out
bampath = args.directory

## functions
### get haplotypes and snp positions from input files
def get_haplotypes(haplotypeFileName):
    
    hapFile = open(haplotypeFileName, 'r')
    haplotypes = {}
    #hList = []
    for line in hapFile:
        line = line.strip('\n')
        parts = line.split('\t')
        haplotypes[parts[1]]=parts[0]
        #hList.append(parts[0])
    hapFile.close()
    
    return haplotypes

def get_variable_pos(positionFileName):
    
    posFile = open(positionFileName, 'r')

    posList = []
    for line in posFile:
        line = line.strip('\n')
        posList.append(int(line))
    posFile.close()
    positions = tuple(posList)
    
    return positions

### get info for each read in bam and counts
def bamInfo(bamFileName,positions,refName):
    bamfileE = pysam.AlignmentFile(bamFileName, "rb")

    hap1 = bamfileE.fetch(refName)

    # collect info for each read
    cigD = {}
    seqD = {}
    for c in hap1:
        qy = c.query_name
        cig = c.cigarstring
        cigRef = cigar_ref(cig)
        # make seqdict with quality
        seq = c.query_sequence
        qual = c.qual
        seqD.setdefault(qy, {})['seq']=seq
        seqD.setdefault(qy, {})['qual']=qual
        # get haplotype for each
        hapStr = ''
        for p in positions:
            posi = posi_finder(cigRef, p)
            if len(seq) > posi:
                base = seq[posi]
                hapStr = hapStr + base
        cigD[qy]=hapStr
    
    return (cigD, seqD)

def hap_counts(cig_dict, haplotype_dict):
    hapCounts = {}
    hapLists = {}
    notFound = []
    for k,v in cig_dict.items():
        if v in haplotype_dict:
            htype = haplotype_dict[v]
            hapLists.setdefault(htype, []).append(k)
        else:
            notFound.append(k)
            hapLists.setdefault('nf', []).append(k)
    for key,value in hapLists.items():
        count = len(value)
        hapCounts[key]=count
    return hapCounts

### output table
def output_table(count_dict, list_haplotypes, tableName='haplotype_count_table.tsv'):
    sample_name_list = list(count_dict.keys())
    tabout = open(tableName, 'w')
    tabout.write('#haplotypes\t')
    for samp in sample_name_list:
        if samp == sample_name_list[-1]:
            tabout.write(samp+'\n')
        else:
            tabout.write(samp+'\t')
    for hap in list_haplotypes:
        tabout.write(hap+'\t')
        for sam in sample_name_list:
            if hap in count_dict[sam]:
                tabout.write(str(count_dict[sam][hap]))
            else:
                tabout.write('0')
            if sam == sample_name_list[-1]:
                tabout.write('\n')
            else:
                tabout.write('\t')
    tabout.close()

## Main function
def main():
    haplo_dict = get_haplotypes(hapFileName)
    positions = get_variable_pos(posFileName)
    reference = refName
    if args.directory:
        path = bampath+'/'
    else:
        path = './'
    # now main loop over bam files
    sampleCounts = {}
    bamlist = open(bamFileList)
    for line in bamlist:
        line = line.strip('\n')
        bamfile = path+'/'+line
        samplename = bamfile.split('/')[-1]
        samplename = samplename.replace('.bam','')
        sampleInfo = bamInfo(bamfile, positions, reference)
        sampleCig = sampleInfo[0]
        sampleCts = hap_counts(sampleCig,haplo_dict)
        sampleCounts[samplename]=sampleCts
    bamlist.close()
    # output to table
    ## get list of haplotypes
    haplotype_list = list(haplo_dict.values())
    haplotype_list.append('nf')
    ## check for table name
    if args.out:
        output_table(sampleCounts,haplotype_list,outfile)
    else:
        output_table(sampleCounts,haplotype_list)

if __name__ == "__main__":
    main()

