#!/usr/bin/env python
from __future__ import division # to be able to do division

#script to calculate synapomorphic characters from sequences from alignment files

from Bio import AlignIO

#then maybe add all alignment files to list
align_file = 'partial_mammal_practice1_aln.fasta'
alignment = AlignIO.read(align_file, 'fasta')

length = alignment.get_alignment_length()

#create a dictionary with sample# : genus name ## note: need to make this part flexible so that different taxonomic levels can be used
name_map = {} # also maybe make a list of each genus, and a dict where key = genus, values = number of species in genus
for i, record in enumerate(alignment):
    fullname = record.id
    newname = fullname.split('_')[0]
    name_map[i]=newname

# use the function to add all alleles to list, then use list length to sort out the monomorphic positions
def number_alleles(column):
    leng = len(column)
    bp_list = []
    for bp in range(0,leng):
        pos = column[bp]
        if pos in bp_list:
            pass
        else:
            bp_list.append(pos)
    #no_alleles = len(bp_list) # maybe take out if a gap
    return bp_list

def informative_alleles(character_dict):
    """a function to filter out uninformative positions"""
    new_dict = {}
    allele_list = []
    for k,v in character_dict.iteritems():
        newd = v
        for key,value in newd.iteritems():
            new_dict.setdefault(key, []).append(value)
    print new_dict
    for keys, values in new_dict.iteritems():
        alleles = set(values)
        print keys
        print alleles
        if len(alleles) > 1:
            allele_list.append(keys)
    return allele_list

pos_list = []
for num in range(0,length):
    position = alignment[:,num] # add num to list
    #print position
    allele_list = number_alleles(position)
    no_alleles = len(allele_list)
    
    if '-' in allele_list:
        no_alleles = no_alleles - 1
    #print no_alleles

    if no_alleles > 1:
        #print position
        pos_list.append(num)
print pos_list
print len(pos_list)

final_chars = {}

for pos in pos_list:
    gen_char_dict = {}
    for sample in range(0,44):
        char = alignment[sample,pos]
        genus = name_map[sample]
        gen_char_dict.setdefault(genus, []).append(char)
    for key, values in gen_char_dict.iteritems():
        #print key
        newlist = gen_char_dict[key]
        check = set(newlist)
        size = len(check)
        if size > 1:
            continue # maybe here add the ambiguous position so that each taxon has a complete set
        else:
            newchar = newlist[0]
            final_chars.setdefault(key, {})[pos]=newchar
        
#print final_chars

# To filter out the characters that have only one alleles; it is variable but only one synapomorphy
## perhaps compare dictionaries in final_chars, if all at one position the same, then delete. 
### try to make it a function 

inform_alleles = informative_alleles(final_chars)
print inform_alleles
print len(inform_alleles)

# the informative_alleles function makes a list of informative alleles, then the next loop creates a filtered dictionary of dicts for each genus
filtered_chars1 = {}
for k,v in final_chars.iteritems():
    tax_dict = v  #final_chars[k]
    #print tax_dict
    for key, value in tax_dict.iteritems():
        if key in inform_alleles:
            filtered_chars1.setdefault(k, {})[key]=value

# now maybe sort the dictionary
print filtered_chars1
print name_map
#now to deal with gaps?, maybe later

# create a genotype for each taxon, that can be searched against



