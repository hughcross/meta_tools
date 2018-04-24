#!/usr/bin/env python

# script to go through sam file, and sort by specific characters

from cigar_funcs import *

samfile = open('sam_hymeno_frax_new.sam', 'r') # note: adjusted gap extension penalty to catch all frax variants
varfile = open('hymen_var1_cts_new.txt', 'w')
stdfile = open('hymen_std_ref_cts_new.txt', 'w')
var1_list = []
var2_list = []
std_ref = []
other_list = []

pos_list = [86,87,113]
geno_dict = {}

for line in samfile:
	if line.startswith('@'):
		continue
	else:
		line = line.strip('\n')
		line_parts = line.split('\t')
		seq_id = line_parts[0]
		cig = line_parts[5]
		if cig.startswith('*'): # in case there are reads with no match
			continue
		else:
			seq = line_parts[9]
			cig_list = cigar_ref(cig)
			for item in pos_list:
				# here have to add exception / if to deal with shorter sequences 
				if len(seq) < item:
					geno_dict.setdefault(seq_id, []).append('-') # adds gap if outside range of seq
				else:
					adjust = pus_finder(cig_list, item)
					base = seq[adjust]
					geno_dict.setdefault(seq_id, []).append(base)

#print geno_dict
# should have dictionary of lists of genotypes, can go by list index to sort now

for key, value in geno_dict.iteritems():
	nulist = geno_dict[key]
	if nulist[0]== 'G':
		if nulist[1]=='C':
			if nulist[2] == 'T':
				var2_list.append(key)
			else:
				std_ref.append(key)
	elif nulist[0] =='C':
		if nulist[1]=='G':
			var1_list.append(key)
	else:
		other_list.append(key)
#print var2_list
print len(var2_list)
#print var1_list
print len(var1_list)
print len(std_ref)
print len(other_list)
var_dict = {}
#petiole2_4july_ion12_mid11_GRUEU:00428:00263
for sample in var1_list:
	sample_parts = sample.split('_')
	new_sample = '_'.join(sample_parts[:2])
	var_dict.setdefault(new_sample, []).append(sample)

#print var_dict.keys()
for k,v in var_dict.iteritems():
	samplist = var_dict[k]
	varfile.write(k+'\t'+str(len(samplist))+'\n')
	#print k
	#print len(samplist)

std_dict = {}
for sample in std_ref:
	sample_parts = sample.split('_')
	new_sample = '_'.join(sample_parts[:2])
	std_dict.setdefault(new_sample, []).append(sample)

#print std_dict.keys()

for k,v in std_dict.iteritems():
	samplist = std_dict[k]
	stdfile.write(k+'\t'+str(len(samplist))+'\n')
	#print k
	#print len(samplist)






