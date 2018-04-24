#!/usr/bin/env python

# sam functions
# note: this version changed to python 3 compatabile on 24/04/18
import re

#note: the following works for usearch mappings that do not use numbers if it is 1
# maybe turn this into a sorted dictionary (has to be sorted)
def cigar_ref(cig):
    if cig and cig[0].isalpha():
        cig = '1'+cig
    cig_mod = re.sub(r'([A-Z])([A-Z])', r'\g<1>1\g<2>', cig)
    ciggy = re.split('[A-Z]', cig_mod) # getting a python3 warning regarding splits, need to fix
    del ciggy[-1]
    ciggyA = re.split('[0-9]*', cig_mod)

    del ciggyA[0]
    #use list comprehension or map function to convert to integers:results = [int(i) for i in results]
    ciggy = [int(i) for i in ciggy]
    return_list = [ciggy, ciggyA]
    return return_list

def cig_ref_degapper(genotype, align_seq): # the genotype does not have to be the full dict, as only need to adjust the keys (positions)
    """to convert genotype from aligned reference seq to account for gaps"""
    gap_pos_list = []
    length = len(align_seq)
    adjusted_geno = {}
    for pos in range(0, length):
        if align_seq[pos] == '-':
            gap_pos_list.append(pos)
    print(gap_pos_list)

    for k,v in genotype.iteritems():
        newpos = k
        print(k)
        for gap in gap_pos_list:
            if k > gap:
                newpos = newpos-1
            else:
                newpos = newpos
        print(newpos)
        adjusted_geno[newpos]=genotype[k]


    return adjusted_geno


def uc_pos_finder(return_list, ref_position): # have to add starting position in
    cig_num = return_list[0]
    cig_alph = return_list[1]
    seq_position = 0
    #seq_pos = start_position - 1
    #if ref_position < seq_pos:
     #   seq_pos = 0
    ## note, not working, have to rethink so that seq_position is additive as it iterates through the 
    
    for code in range(0,len(cig_alph)):
        if cig_alph[code] == 'M':
            if ref_position < int(cig_num[code]) + seq_position:
                if code == 0:
                    seq_position = ref_position
                else:
                    seq_position = ref_position - seq_position 
        elif cig_alph[code] == 'I':
            seq_position = seq_position - int(cig_num[code])
        elif cig_alph[code] == 'D':
            seq_position = seq_position + int(cig_num[code])
    return seq_position

def pos_finder(return_list, ref_position): # have to add starting position in
    cig_num = return_list[0]
    cig_alph = return_list[1]
    seq_position = 0
    #seq_pos = start_position - 1
    #if ref_position < seq_pos:
     #   seq_pos = 0
    ## note, not working, have to rethink so that seq_position is additive as it iterates through the 
    
    for code in range(0,len(cig_alph)):
        if cig_alph[code] == 'M':
            if ref_position <= cig_num[code] + seq_position:
                if code == 0:
                    seq_position = ref_position
                else:
                    seq_position = ref_position - seq_position 
        elif cig_alph[code] == 'I':
            seq_position = seq_position + int(cig_num[code])
        elif cig_alph[code] == 'D':
            seq_position = seq_position - int(cig_num[code])
    seq_position = seq_position -1
    return seq_position

def pus_finder(return_list, ref_position): # have to add starting position in
    cig_num = return_list[0]
    cig_alph = return_list[1]
    for num in range(1,len(cig_num)+1):
        tot = sum(cig_num[:num])
        if ref_position <= tot:
            listnum = num - 1
            break
        else:
            listnum = num # changed from pass
    for index in range(0,listnum):
        if cig_alph[index] == 'M':
            ref_position == ref_position
        elif cig_alph[index]== 'I':
            ref_position = ref_position + cig_num[index]
        elif cig_alph[index] == 'D':
            ref_position = ref_position - cig_num[index]
        elif cig_alph[index] == 'S':
            ref_position = ref_position + cig_num[index]
    adj_pos = ref_position - 1

    return adj_pos












