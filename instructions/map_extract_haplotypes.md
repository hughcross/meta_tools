# How to use the *extractHaplotypes* script

This page goes through the steps to prepare the data to extract the haplotype information

## Prepping the data

After you have trimmed your raw reads for quality and removed the adapters (using AdapterRemoval or Cutadapt), you have to merge the paired ends into single reads (you can run this in AdapterRemoval as part of the trimming) and then trim the ends so it can align well to the primer (if you haven't done this with AR or cutadapt). You can use a tool like Pear or BBmerge to do this. Here is an example using BBmerge. 

On NeSI load the BBMap and seqtk modules:

```
module load BBMap
module load seqtk/1.3-gimkl-2018b
```

Then you can run these commands on all your samples in a loop:

```
cat samplelist | while read samp;
do
bbmerge.sh \
 in1=${samp}.pair1.truncated \
 in2=${samp}.pair2.truncated \
 out=${samp}.merged.fastq

seqtk trimfq -b 22 -e 20 ${samp}.merged.fastq > ${samp}.trimd.fastq

done
```

## Map reads to haplotype sequence

Next you have to map all these reads to a haplotype sequence. I used haplotype 1 but it doesn't matter. Most reads will map and to this sequence and then record where it varies from the reference sequence. 

First, index the reference sequence:

```
module load BWA

bwa index haplotype_ref.fasta
```

Then you want to map each sample to the reference, then sort and index the map file. You can do this all in a loop. I suggest making a new folder for the mapping files, and then use a variable for the reference and data folders:

```
DATADIR='/path/to/data/folder/'
REFDIR='path/to/reference/sequence/'
```

Then run the loop (make sure a list of samples is in the mapping folder; you can just copy it from the data folder from the above step)

```
cat samplelist | while read samp;
do
bwa mem -t 4 ${REFDIR}/synthHaplotype1.fasta ${DATADIR}/${samp}.trimd.fastq > ${samp}.sam
samtools view -b -@ 4 ${samp}.sam | samtools sort -@ 4 - -o ${samp}.bam
samtools index ${samp}.bam
done
```

## Run the Python script to get counts of all haplotypes

To run the script, you just need to include two text files: one that lists the positions of each SNP, and another, tab-delimited file that lists the haplotype name, followed by the bases at each of the variable positions. There are examples in this repository of these files. Once you have these, you can proceed with the script.

Here shown as a loop through all the sample bam files

```
cat sample2list | while read samp;
do
/path/to/github/repos/meta_tools/extractHaplotypes.py \
  -b ${samp}.bam \
  -p /path/to/positions/file/june2020_positions.txt \
  -H /path/to/haplotype/file/paua_bucket_haps.txt \
  -r Haplotype1 \
  -o haplotypes_${samp}
done
```

The `-r` option in the script is the sequence name of the reference sequence to which each sample sequence file was mapped. This is not the name of the reference file, but the name of the sequence itself (the text after the > and before the first space). 

The above script will create a folder for each sample that contains fasta files and a logfile with the counts for each haplotype. You could copy all the logfiles into one folder for easy access:

```
mkdir logfiles

cat samplelist | while read samp;
do
cp haplotypes_${samp}/${samp}.log logfiles
done
```

(Note: above assumes that each sample folder is named as in the previous script, e.g. haplotypes_sample1)

## Checking for noise in the haplotype output

As the Python script is only looking at the loci that define the haplotypes (in the present example, seven). Therefore, it might be good to check on what kind of variability there is at the other positions. You can do this for several sample haplotype files, or all if you like. Here is the code for a single sample (in this example the reads for haplotype 1 from Sample-10_7_1):

```
module load VSEARCH/2.4.3-gimkl-2017a

vsearch --derep_fulllength Sample-10_7_1_haplotype1.fasta \
  --minuniquesize 10 --relabel rep. --sizeout \
  --output rl_dr10_Sample-10_7_1_haplotype1.fasta

```

You can open this output in Geneious to see the variation. The 'size=' after each sequence name will tell you how many of each replicate there is. In Geneious, right click on the sequence view and select sort > by name, so that the biggest replicates are at the beginning. 








