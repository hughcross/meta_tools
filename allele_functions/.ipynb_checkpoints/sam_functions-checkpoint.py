#!/usr/bin/env python

# functions to manipulate sam files
from __future__ import division
# make a function that opens a single sam file, and outputs a dictionary of sequences matching list of contigs
## dict: {contig1 : [read1, read2, read3], contig2 : [read4, read5, read6]}
### maybe better have dictionary of dictionaries: {contig1 : {readname1 : read1, readname2, read2}, contig2 : {} ...}
import subprocess
import sys
import select
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import io #from Bio.SeqRecord import SeqRecord
# the following changed for python 3 (24/04/18)
#from io import io #from StringIO import StringIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from cigar_funcs import cigar_ref
from collections import OrderedDict


def grab_sam_reads(samfile, contiglist):
    """function to grab reads mapping to all contigs in a list, output a dict of dicts"""
    read_dict = {}
    for line in samfile:
        if line.startswith('@'):
            continue
        else:
            line = line.strip('\n')
            line_parts = line.split('\t')
            read = line_parts[0]
            contig = line_parts[2]
            cig = line_parts[5]
            readMD = read+';'+cig ## meaning 'read with MetaData'
            if contig in contiglist:
                #maq = line_parts[4]
                seq = line_parts[9]
                #read_dict.setdefault(contig, []).append(seq) # create a list of seqs matching that contig
                read_dict.setdefault(contig, {})[readMD]=seq
                ## maybe add metadata later
                #if maq > mq:
                #    if len(seq) > minseq:
                #        seqdict.setdefault(read, seq)
                #        condict.setdefault(contig, []).append(read)
    return read_dict

## need to convert dictionary of reads to biopython sequence object
def seq_dict_to_seqIO(dict_name):
    """convert a sequence dictionary to biopython sequence objects"""
    handle = io.StringIO()
    for key, value in dict_name.iteritems():
        rec = SeqRecord(Seq(value, IUPAC.unambiguous_dna), id=key, description='') # maybe add metadata to description later
        SeqIO.write(rec, handle, 'fasta') # this writes each sequence to a variable holding data
    data = handle.getvalue()
    return data

def seq_dict_to_seqIO_no_overhangs(dict_name):
    """convert a sequence dictionary to biopython sequence objects, filtering out overhangs from sam file mapping"""
    handle = io.StringIO()
    for key, value in dict_name.iteritems():
        full_seq = value
        full_id = key
        read_name = full_id.split(';')[0]
        cigar_score = full_id.split(';')[1]
        full_cigar = cigar_ref(cigar_score)
        cigN = full_cigar[0]
        cigT = full_cigar[1]
        if cigT[0] == 'S':
            new_seq = full_seq[cigN[0]:]
        else: 
            new_seq = full_seq
        if cigT[-1] == 'S':
            newer_seq = new_seq[:-cigN[-1]]
        else:
            newer_seq = new_seq
        # maybe modify later to include metadata in output under description to indicate if it has been trimmed
        rec = SeqRecord(Seq(newer_seq, IUPAC.unambiguous_dna), id=read_name, description='') # maybe add metadata to description later
        SeqIO.write(rec, handle, 'fasta') # this writes each sequence to a variable holding data
    data = handle.getvalue()
    return data

def seq_dict_rm_overhangs(dict_name):
    """remove cig overhangs from sequence dictionary and output new sequence dictionary"""
    #handle = StringIO()
    rm_dict = {}
    for key, value in dict_name.iteritems():
        full_seq = value
        full_id = key
        read_name = full_id.split(';')[0]
        cigar_score = full_id.split(';')[1]
        full_cigar = cigar_ref(cigar_score)
        cigN = full_cigar[0]
        cigT = full_cigar[1]
        if cigT[0] == 'S':
            new_seq = full_seq[cigN[0]:]
        else: 
            new_seq = full_seq
        if cigT[-1] == 'S':
            newer_seq = new_seq[:-cigN[-1]]
        else:
            newer_seq = new_seq
        rm_dict[key]=newer_seq
        # maybe modify later to include metadata in output under description to indicate if it has been trimmed
        #rec = SeqRecord(Seq(newer_seq, IUPAC.unambiguous_dna), id=read_name, description='') # maybe add metadata to description later
        #SeqIO.write(rec, handle, 'fasta') # this writes each sequence to a variable holding data
    #data = handle.getvalue()
    return rm_dict

def get_muscle_align(data): 
    muscle_cline = MuscleCommandline() # input=data, 
    stdout, stderr = muscle_cline(stdin=data)
    alignment = AlignIO.read(io.StringIO(stdout), 'fasta')
    # temporary output to file
    #print(alignment)
    return alignment

def get_align_rep_freq(align):
    """to take derepped alignment and output list of frequencies for each rep"""
    freq_list = []
    for record in align:
        id = record.id
        size = id.split(';')[-2][5:]
        freq_list.append(int(size))
    return freq_list

# format for reads in alignment
##Ac_argyrophylla_85_P3Z9O:00659:00693;45S64M1I21M2I11M;size=1;
def OLDget_haplotypes(alignment,list_of_var_pos):
    """to obtain haplotypes from alignment by position with freq in name"""
    final_gt_list = []
    allele_dict = OrderedDict()
    num_reads = len(alignment)
    newstr = ''
    count = 0
    gt_dict = {}
    for i in range(0,num_reads): # number of sequences ## change this to parsing the Alignment
        name = 'gt'+str(count)
        # id = record.id # then parse this to get depth back out
        # seq = record.seq
        for pos in list_of_var_pos:
            geno = alignment[i,pos]
            # add an if to remove gaps from alleles
            gt_dict.setdefault(name, []).append(geno) # maybe make this a dictionary of dicts, to include haplotype depth
            count += 1

    for k,v in gt_dict.iteritems():
        if v in final_gt_list:
            continue
        elif '-' in v: # try to prevent gaps being included in alleles ### may need to change this for variable length algnments
            continue
        else:
            final_gt_list.append(v)
    nu_count = 1

    for al in final_gt_list:
        nu_name = 'hap'+str(nu_count)
        values = gt_dict.values() #####
        freq = values.count(al) # need to change this to get right number of each haplotype
        freq_name = nu_name+'_freq='
        allele_dict[nu_name]=al
        allele_dict[freq_name]=freq
        nu_count += 1
    hap_numplus = len(final_gt_list)+1
    #newstr = ''
    for h in range(1,hap_numplus):
        typ = 'hap'+str(h)
        freeq = typ+'_freq='
        current = typ+'='+str(allele_dict[typ])+':freq='+str(allele_dict[freeq]) 
        if h == 1:
            newstr = newstr+current
        else:
            newstr = newstr+'; '+current
    num_alleles = len(final_gt_list)
    output = (num_alleles,newstr)
    return output

def get_haplotypes(alignment,list_of_var_pos):
    """to obtain haplotypes from alignment by position with freq in name"""
    final_gt_list = []
    allele_dict = OrderedDict()
    num_reads = len(alignment)
    newstr = ''
    count = 0
    gt_dict = {}
    rep_freq = {}
    for record in alignment:
        seq_id = record.id
        #size = seq_id.split(';')[-2][5:]
        read = record.seq
        #name = 'gt'+str(count)
        gt_str = ''
        # id = record.id # then parse this to get depth back out
        # seq = record.seq
        for pos in list_of_var_pos:
            geno = read[pos]
            gt_str = gt_str+geno

        gt_dict[seq_id]=gt_str
        count += 1

    final_gt_dict = {}
    for k,v in gt_dict.iteritems():
        size = int(k.split(';')[-2][5:])
        if '-' in v:
            continue
        else:
            final_gt_dict.setdefault(v, []).append(size)
    #print final_gt_dict

    nu_count = 1
    for key, value in final_gt_dict.iteritems():
        
        freq = sum(value)
        if freq > 0: # changed to check on these
            nu_name = 'hap'+str(nu_count)
            freq_name = nu_name+'_freq='
            allele_dict[nu_name]=key
            allele_dict[freq_name]=freq
            nu_count += 1
    hap_numplus = len(allele_dict)+1
    #newstr = ''
    for h in range(1,hap_numplus):
        typ = 'hap'+str(h)
        freeq = typ+'_freq='
        if typ in allele_dict:
            current = typ+'='+str(allele_dict[typ])+':freq='+str(allele_dict[freeq]) 
            if h == 1:
                newstr = newstr+current
            else:
                newstr = newstr+';'+current
    num_alleles = int(len(allele_dict)/2)
    output = (num_alleles,newstr)
    return output



