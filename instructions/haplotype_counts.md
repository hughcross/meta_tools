# Creating table of counts for multiple samples

These instructions are to get a single table from a list of bam files. This will not produce fasta files for each sample, as in the `extractHaplotypes.py` script. 

## Prepping the samples

To prep the samples, follow the same instructions as for the [**extract haplotypes instructions**](map_extract_haplotypes.md), including **Prepping the data** and **Map reads to haplotype sequence**.

The only other thing you need is a file with a list of the bam files to be used to make the table, one per line. 

## Run the script


```
module load Python/3.8.2-gimkl-2020a

SCRIPTDIR='/nesi/project/uoo02328/programs/meta_tools/scripts'

# to get help
${SCRIPTDIR}/haplotype_counts.py -h

# example 1
${SCRIPTDIR}/haplotype_counts.py \
 -b bamlist \
 -p dr100_positions.txt \
 -H dr100__haplotypes.txt \
 -r Haplotype1

# example 2: bam files in another folder, naming output
${SCRIPTDIR}/haplotype_counts.py \
 -b bamlist \
 -p dr100_positions.txt \
 -H dr100__haplotypes.txt \
 -r Haplotype1 \
 -d /path/to/bam/files \
 -o table_output_name
```

Where `-p` and `-H` are the outputs from the `find_haplotypes.py` script, `-b` is the file listing the bam files, `-r` is the name of the reference sequence (the sequence name in the fasta file), `-d` is the path to the folder with bam files, and `-o` is a name for the output table. 