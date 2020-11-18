#!/usr/bin/env python


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
    for k,v in character_dict.items():
        newd = v
        for key,value in newd.items():
            new_dict.setdefault(key, []).append(value)
    #print new_dict
    for keys, values in new_dict.items():
        alleles = set(values)
        #print keys
        #print alleles
        if len(alleles) > 1:
            allele_list.append(keys)
    return allele_list

