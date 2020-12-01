#!/usr/bin/env python
# coding: utf-8

# # Finding haplotypes from alignment

from Bio import AlignIO
from Bio import SeqIO

## cmd line args
import sys
import argparse
# arguments

parser = argparse.ArgumentParser(description='find haplotypes in alignment files')

parser.add_argument('-f', '--aln_file', dest='aln_input',
type=str,
help='alignment input in fasta format')

parser.add_argument('-p', '--prefix', dest='output_prefix',
type=str,
help='prefix for all output files [optional]')

parser.add_argument('-m', '--minor_allele_freq', dest='minor_allele_frequency',
type=str,
help='minimum minor allele frequency to call SNP [optional: default=0.1]')

parser.add_argument('-H', '--hap_freq', dest='minor_hap_freq',
type=str,
help='minimum minor haplotype frequency to include haplotype [optional: default=0.1]')

### argparse parameters
args = parser.parse_args()
align = args.aln_input
outprefix = args.output_prefix
maf = args.minor_allele_frequency
mhf = args.minor_hap_freq

## set defaults for options
if args.output_prefix:
    prefix = outprefix
else:
    prefix = 'found'

if args.minor_allele_frequency:
    mmaf = maf
else:
    mmaf = 0.1

if args.minor_hap_freq:
    mmhf = mhf
else:
    mmhf = 0.1

# ## make functions
def global_adjust(alignment):
    """this assumes that the reference is the first sequence"""
    alnLength = alignment.get_alignment_length()
    for p in range(0,alnLength):
        if alignment[0,p] == '-':
            continue
        else:
            startpoint = p
            break
    for e in range(1,alnLength):
        if alignment[0,-e] == '-':
            continue
        else:
            endpoint = e-1
            break
    new_alignment = alignment[:,startpoint:-endpoint]
    return new_alignment

def sample_info(alignment):
    sampleDict = {}
    sampleCounts = {}
    sampleNames = []
    for samp in range(1,len(alignment)):
        sequence = str(alignment[samp,:].seq).upper()
        seqID = alignment[samp,:].id
        sampleNames.append(seqID)
        counts = int(seqID.split('size=')[1])
        sampleDict[seqID]=sequence
        sampleCounts[seqID]=counts
    return (sampleDict, sampleCounts, sampleNames)

def pos_counts(sampleDict, sampleCounts, sampleNames, alnLength):
    posCounts = {}
    for pos in range(0,alnLength):
        baseDict = {}
        for sample in sampleNames:
            base = sampleDict[sample][pos]
            ct = sampleCounts[sample]
            baseDict.setdefault(base,[]).append(ct)
        totalCounts = {}
        for k,v in baseDict.items():
            total = sum(v)
            totalCounts[k]=total
        posCounts[pos]=totalCounts
    return posCounts

def informative_alleles(posCounts, mmaf=0.1):
    infAlleles = []
    for k,v in posCounts.items():
        ctlist = list(v.values())
        topCount = max(ctlist)
        sumCounts = sum(ctlist)
        for key,value in v.items():
            if value != topCount:
                if (value/sumCounts) > mmaf:
                    infAlleles.append(k)
    return infAlleles

def informative_bases(infAlleles, posCounts):
    infBases = {}
    for i in infAlleles:
        baseList = list(posCounts[i].keys())
        bases = set(baseList)
        infBases[i]=bases
    return infBases

def raw_haplotypes(sampleDict, sampleCounts, infAlleles, infBases):
    haplotypes = {} # hap:count
    for k,v in sampleDict.items():
        hapStr = ''
        for i in infAlleles:
            base = v[i]
            if base in infBases[i]:
                hapStr += base
            else:
                hapStr += '-'
        haplotypes.setdefault(hapStr, []).append(sampleCounts[k])
    return haplotypes

def major_haplotypes(haplotypes, mhf=0.1):
    """mhf == major haplotype frequency"""
    hapTotals = {}
    for k,v in haplotypes.items():
        total = sum(v)
        hapTotals[k]=total
    totalFreq = sum(hapTotals.values())
    major_haps = []
    for key,value in hapTotals.items():
        if value/totalFreq >= mhf:
            major_haps.append(key)
    return major_haps

def output_files(prefix, major_haplotypes, informative_alleles):
    positionFileName = prefix+'_positions.txt'
    haplotypeFileName = prefix+'_haplotypes.txt'
    posout = open(positionFileName, 'w')
    hapout = open(haplotypeFileName, 'w')
    counter = 1
    last_hap = major_haplotypes[-1]
    last_pos = informative_alleles[-1]
    for hap in major_haplotypes:
        if hap == last_hap:
            hapout.write(str(counter)+'\t'+hap)
        else:
            hapout.write(str(counter)+'\t'+hap+'\n')
            counter += 1
    for p in informative_alleles:
        if p == last_pos:
            posout.write(str(p+1))
        else:
            posout.write(str(p+1)+'\n')
    posout.close()
    hapout.close()

## added to remove gaps in reference
def remove_ref_gaps(positions, alignment):
    """remove gaps from reference to adjust informative positions"""
    ref = str(alignment[0,:].seq)
    adjA = []
    for pos in positions:
        sub = ref[:pos]
        nogap = sub.replace('-','')
        #print(pos, len(sub), 'revised:',len(nogap)+1)
        adjust = len(nogap)
        adjA.append(adjust)
    return adjA

## Main function
def main():
    alnFile = align
    aln1 = AlignIO.read(alnFile,'fasta')

    #aln1.get_alignment_length()
    aln2 = global_adjust(aln1)
    alnLength = aln2.get_alignment_length()

    aln2_info = sample_info(aln2)
    sampleNames = aln2_info[2]
    aln2_pos_counts = pos_counts(aln2_info[0], aln2_info[1], aln2_info[2], alnLength)
    aln2_infA = informative_alleles(aln2_pos_counts)
    aln2_infB = informative_bases(aln2_infA, aln2_pos_counts)
    aln2_rh = raw_haplotypes(aln2_info[0], aln2_info[1], aln2_infA, aln2_infB)
    aln2_mh = major_haplotypes(aln2_rh)
    aln2_adj = remove_ref_gaps(aln2_infA, aln2)

    output_files(prefix, aln2_mh, aln2_adj)

if __name__ == "__main__":
    main()

