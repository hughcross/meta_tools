# Find haplotypes from amplicon sequence data

The following assumes you have single or paired end reads amplicon sequence reads that have been quality-filtered, merged (for PE reads), and primer sequence has been clipped off (although should still work if primer sequence is there). Also it will help to set minimum and maximum lengths for your output reads (the min. should be only about 10 - 20 bp less than the reference sequence).

## Combine all reads

The first thing is to combine all merged, filtered reads from all samples into one file (similar to what you do to find OTUs).

```
DATA='/path/to/data'

cat ${DATA}/*.trimd.fastq > allSamples.fastq
```

Convert the fastq to a fasta file:

```
module load seqtk/1.3-gimkl-2018b

seqtk seq -A allSamples.fastq > allSamples.fasta
```

## Dereplicate the sequence file

Now dereplicate the file with all reads. 


```
module load VSEARCH/2.14.2-GCC-9.2.0

vsearch --derep_fulllength allSamples.fasta \
  --minuniquesize 100 --relabel rep. --sizeout \
  --output dr100_allSamples.fasta
```

Check the number of uniques written. The final readout will say something like:

*2289 uniques written, 267566 clusters discarded (99.2%)*

If the number of uniques is more than three or four thousand, set the `--minuniquesize` parameter higher. Try 1000 or more. 

### Add the reference sequence to the dereplicate output fasta

This is the reference sequence that you will map the reads to. The reference fasta should contain only a single sequence.

```
cat haplotype_ref.fasta dr100_allSamples.fasta > dr100_allSamples_wRef.fasta
```

Note: make sure that the reference fasta is first in above command, so it is the first sequence in the resulting file

## Align reference and unique sequences

Now align the sequences. This can take a few minutes, depending on the number of sequences and length (see note under vsearch command).

```
module load MAFFT/7.429-gimkl-2020a

mafft --auto dr100_allSamples_wRef.fasta > aln_dr100_allSamples_wRefa.fasta
```

## Run Python script to find haplotypes

Now, run the script

```
module load Python/3.8.2-gimkl-2020a

SCRIPTDIR='/nesi/project/uoo02328/programs/meta_tools/scripts'

# for help
${SCRIPTDIR}/find_haplotypes.py -h 

# example run
${SCRIPTDIR}/find_haplotypes.py -f aln_dr100_allSamples_wRef.fasta -p dr100
```

The output of this script will produce two files: a file indicating haplotypes and another with positions of SNPs. These are the two files needed for either the `extractHaplotypes.py` or `haplotype_counts.py` scripts.


