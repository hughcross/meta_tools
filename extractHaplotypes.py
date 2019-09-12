#!/usr/bin/env python

# coding: utf-8
## todo: create sequence files from unidentified reads (no haplotype found)
# # Extract Haplotype sequences from bam file
import sys
import os
import pysam
from allele_functions.cigar_funcs import *
import argparse
# ## Get variables needed (use Argparse)

parser = argparse.ArgumentParser(description='extract haplotypes from bam files')

parser.add_argument('-b', '--bam_file', dest='input',
type=str,
help='bam input filename')

parser.add_argument('-p', '--positions', dest='posifile',
type=str,
help='file with haplotype SNP positions on alignment')

parser.add_argument('-H', '--haplotypes', dest='haplos',
type=str,
help='tab-delimited file with haplotype name (1st col) and hap string (2nd col)')

parser.add_argument('-r', '--ref_name', dest='ref',
type=str,
help='Reference Name')

parser.add_argument('-f', '--format', dest='format',
type=str,
help='output format, either fastq or fasta (default)')

parser.add_argument('-o', '--output_folder', dest='out',
type=str,
help='output folder name (default: name of sample')

### argparse parameters
args = parser.parse_args()
bamFileName = args.input
hapFileName = args.haplos
posFileName = args.posifile
refName = args.ref
outfolder = args.out
outformat = args.format

####
#hapFileName = '/Users/hughcross/Analysis/repos/meta_tools/example_haplotypes.txt'
#posFileName = '/Users/hughcross/Analysis/repos/meta_tools/example_position_list.txt'
#bamFileName = '/Users/hughcross/Analysis/Paua/sort_allBucks.bam'
#refName = 'CCB11_Haplotype1'

bamNameNoPath = bamFileName.split('/')[-1]
sampleName = bamNameNoPath.split('.')[0]
print('sample name:', sampleName)

logFileName = sampleName+'.log'
print('log file:', logFileName)

#### make the output dir
if args.out:
    os.makedirs(outfolder)

## logfile
if args.out:
    logFileName = outfolder+'/'+sampleName+'.log'
else:
    logFileName = sampleName+'.log'

# ## Get dict of haplotypes and list of positions from files
hapFile = open(hapFileName, 'r')

hapDict = {}
hList = []
for line in hapFile:
    line = line.strip('\n')
    parts = line.split('\t')
    hapDict[parts[1]]=parts[0]
    hList.append(parts[0])
hapFile.close()
print('number of haplotypes:', len(hapDict))
print('haplotypes', hapDict)
#print(hList)

posFile = open(posFileName, 'r')

posList = []
for line in posFile:
    line = line.strip('\n')
    posList.append(int(line))
posFile.close()
positions = tuple(posList)
#print(posList)
print('SNP positions in alignment', positions)

# ## Use pysam to parse the file

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
#print(len(cigD))
#print(len(seqD))


hapLists = {}
notFound = []
for k,v in cigD.items():
    if v in hapDict:
        htype = hapDict[v]
        hapLists.setdefault(htype, []).append(k)
    else:
        notFound.append(k)
        
print('\nNumber of unidentified reads:', len(notFound))

# ## Output to sequence files

# for fasta file

# for fastq file
if outformat == 'fastq':
    for hap in hList:
        if args.out:
            newFileName = outfolder+'/'+sampleName+'_haplotype'+hap+'.fastq'
        else:
            newFileName = sampleName+'_haplotype'+hap+'.fastq'
        print('writing:', newFileName)
        output = open(newFileName, 'w')
        readlist = hapLists[hap]
        for read in readlist:
            output.write('@'+read+'\n'+seqD[read]['seq']+'\n+\n'+str(seqD[read]['qual'])+'\n')
        output.close()
else: # makes fasta default
    for hap in hList:
        if args.out:
            newFileName = outfolder+'/'+sampleName+'_haplotype'+hap+'.fasta'
        else:
            newFileName = sampleName+'_haplotype'+hap+'.fasta'
        print('writing:', newFileName)
        output = open(newFileName, 'w')
        readlist = hapLists[hap]
        for read in readlist:
            output.write('>'+read+'\n'+seqD[read]['seq']+'\n')
        output.close()

# ## Output logfile
print('writing logfile', logFileName)
logout = open(logFileName,'w')

logout.write('Haplotype\tNumber of reads\n')

for h in hList:
    numReads = str(len(hapLists[h]))
    logout.write(h+'\t'+numReads+'\n')
logout.write('None\t'+str(len(notFound))+'\n')
logout.close()

