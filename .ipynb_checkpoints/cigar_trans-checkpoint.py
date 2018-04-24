#!/usr/bin/env python

# script to identify positions with cigar scores

from allele_functions.cigar_funcs import *

cig_score = raw_input('enter cigar score: ')

pos = int(raw_input('enter position to check: '))

seq = raw_input('enter the sequence: ')

cig_list = cigar_ref(cig_score)


print cig_list[0]
print cig_list[1]


adjustment = pus_finder(cig_list, pos)

print adjustment


# so far, now have to apply adjustment to sequence (can input for now), then spit out base at that position

#final_pos = adjustment + pos

base = seq[adjustment]

print base



