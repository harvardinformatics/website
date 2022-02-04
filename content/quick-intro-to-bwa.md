Title: Short Introduction to bwa
Date: 2020-04-06
Category: Tutorials
Author: Nathan Weeks
Tags: grep, unix
Summary: This tutorial provides a basic overview of using bwa to align fastqs to a reference genome

[TOC]

## Aligning Reads to a Reference Genome with BWA
### Bioinformatics Coffee Chat - April 7, 2020
### Programs used in this session
#### [BWA](http://bio-bwa.sourceforge.net/), [Picard](https://broadinstitute.github.io/picard/)

##### You can run this example yourself [here](https://github.com/harvardinformatics/intro-to-bwa)

Let's get started!

We have two directories:
* `00_genome` - reference sequence
* `01_fastqs` - raw Illumina reads to map to the genome

#### Fastq -> SAM -> BAM

#### Process reference genome
BWA requires building an index for your reference genome to allow it to more efficiently search the genome during sequence alignment:


```bash
bwa index -p 00_genome/Falb 00_genome/Falbicolis.chr5.fa.gz
```

You should have several new files in the `00_genome` directory that all start with 'Falb', since this is the value we gave after the `-p` flag.

Lets use BWA to align our reads to the reference genome and get a SAM file (Sequence Alignment/Map), one for each sample. For information about BWA, simply type its name, followed by `--help`. This procedure works for most programs.


```bash
mkdir -p 02_bams
for INDEX in 1 2 3 34 35;
do
   bwa mem -M -t 1 -R "@RG\tID:COL_${INDEX}\tSM:COL_${INDEX}" 00_genome/Falb \
   01_fastqs/Falb_COL${INDEX}.1.fastq.gz \
   01_fastqs/Falb_COL${INDEX}.2.fastq.gz \
   > 02_bams/Falb_COL${INDEX}.sam #2> Falb_COL${INDEX}.log
done
```

Note: read groups refer to sets of reads that were generated from a single run of a sequencing instrument, and this information is stored in our SAM file on lines that start with `@RG`.
Using read groups allows us to not just distinguish between samples, but also particular samples that were sequenced across several experiments.
Programs like the GATK require this information so that it can attempt to compensate for variability between sequencing runs.

Next, lets convert our SAM files from BWA to BAM files, which are compressed versions that a lot of downstream programs use as input files.


```bash
for INDEX in 1 2 3 34 35;
do
  picard SortSam \
  I=02_bams/Falb_COL${INDEX}.sam \
  O=02_bams/Falb_COL${INDEX}.sorted.bam \
  SORT_ORDER=coordinate \
  CREATE_INDEX=true
done
```

Last, we need to create an index of our BAM file in order for downstream programs to quickly access its contents.


```bash
for INDEX in 1 2 3 34 35
do
  picard BuildBamIndex \
  I=02_bams/Falb_COL${INDEX}.sorted.bam
done
```

### That's all, folks!
