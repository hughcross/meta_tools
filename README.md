# meta_tools

A set of tools for working with metagenomic data. 

## Population genomics of fur seals 

The scripts for this project are used to determine the haplotypes of the mitogenome assemblies from an environmental bait capture experiment of New Zealand fur seals. 

### Overall strategy 

Because the assemblies are derived from environmental DNA, it is not possible to assign the entire assembly to a single haplotype, as multiple individuals could contribute to any single water sample. Therefore, the strategy I used was to assign **each sequence read** to a haplotype, as best as possible. 

The following alignment from part of the **cytB** mitochondrial gene illustrates the approach:

![fur seal alignment](images/furseal_haplo_finder_slide1_lg.png)


The alignment shows the consensus alignment of 20 **cytB** fur seal haplotypes, with reads (A) - (H) aligned below. The included scripts will assign each read to a haplotype (or haplotypes) depending on their position:

read (A) is assigned to haplotype 15

read (B) is assigned to haplotype 14

read (C) is assigned to haplotype 19

read (D) is assigned to haplotype 4

The remaining reads are less certain:

read (E) could be haplotype 8, 9, 10, 13, 17, 18, or 20

read (F) looks to be haplotype 19, but has a SNP (at position 533) that is not in any of the known haplotypes

read (G) could be haplotypes 1, 2, 3, 5, 6, 7, 11, or 12

read (H) has many SNPs not found in any of the haplotypes, and is probably sea lion or another species.

The reads aligned to the haplotypes in this example are made up, just to illustrate the different situations. From this example, reads A - D can be assigned clearly. Reads E and G could be one of several haplotypes so are indeterminate. Read F is probably haplotype 19, though it has an additional SNP. Whether or not to count this SNP depends on its frequency. If the SNP on position 533 has a very low frequency, it can probably be discounted. If it has a high frequency, it may be a new variant of haplotype 19, however, as the script was set up to count frequencies of existing haplotypes only, it is an open question how to deal with these. We could just delete any read that is not a perfect match to any known haplotypes. 

Read H looks to be completely different. Something that far off any of the known haplotypes should probably discounted. Reads that are very far from any of the haplotypes would likely not map to the mitogenome at all. This is a good reason to check what is aligning to sea lion (or other relatives).

The scripts as they are just count the number of reads associated with any haplotype. This could be made more rigorous by assigning a confidence to each read by number of informative sites, or similar. 

Suggestions welcome! 



----------------------

This covers both amplicon/single gene metagenomics (metabarcoding), and shotgun/functional metagenomics

[**Instructions for running *find_haplotypes.py* script**](instructions/finding_haplotypes.md)

[**Instructions for running *haplotype_counts.py* script**](instructions/haplotype_counts.md)

[**Instructions for running *extractHaplotypes.py* script**](instructions/map_extract_haplotypes.md)

