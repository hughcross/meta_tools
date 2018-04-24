hymenoscyphus remapping


first, grep out from qiime taxonomy file, all Hymenoscyphus matches,

then short python script to parse usearch map file to id all reads matching Hymen. then
extract these reads from ngs file to new file.

	usearch_map_parser.py
	
then use ipython to extract reference seqs from main ref fasta. 
from these files:
hymenoscyph_refs.txt
extracted_new_ion_refs_unite_ncbi_its.fasta
and from
extracted_ion_ash_swarms_full.fasta

combine:
cat hymenoscyph_ion_refs.fasta hymenoscyph_ion_swarms.fasta > hymenoscyph_ion_its_seqs.fasta

rename sequences with script:

newfile
hymenoscyph_ion_its.fasta

new_combined_ion_tax_qiime_format.txt


import to Geneious, make alignment, 


use alignment to design script to assign reads to H. spp and H. frax variants (using Cigar scores)

now try basic to distinguish among H frax vars
first map using bwa
bwa index -p Hymenoscyphus_fraxineus_varA_58_its2 Hymenoscyphus_fraxineus_varA_58_its2.fasta 
bwa mem -E 2[4] Hymenoscyphus_fraxineus_varA_58_its2 hymenoscyphus_read_matches.fasta > sam_hymeno_frax_new.sam

178 'T' variant (variant 2)
2073 'CG' variant (variant 1)
6798 'GC' variant (standard variant)
20 others
26 did not match anything

redo tables in excel
