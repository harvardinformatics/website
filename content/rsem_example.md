Title: RSEM example on Odyssey
Author: Adam Freedman
Date: 2016-01-04 00:00
Category: Software
Tags: Next-Gen Sequencing, Transcriptome, RNA-seq Quantitation, Differential Expression, RSEM
Summary: An example of quantifying RNA-seq expression with RSEM on Odyssey cluster


[RSEM](http://deweylab.github.io/RSEM/README.html) is a software package for estimating gene and isoform expression levels from single-end or paired-end RNA-Seq data. The software works with transcriptome sequences and does not require a reference genome. It can either perform the read alignment step prior to quantification, or take an alignment (bam) file as input, so long as the alignment settings are appropriate for RSEM. Currently, RSEM can perform the alignment step with three different aligners: bowtie, bowtie2, or STAR. It uses the Expectation Maximization (EM) algorithm to estimate abundances at the isoform and gene levels.

If an annotated reference genome is available, RSEM can use a gtf file representation of those annotations to extract the transcript sequences for which quantification will be performed, and build the relevant genome and transcriptome indices. If this is done and the alignment step is implemented within RSEM, the option is available to also write the read alignments in genomic coordiinates, permitting visualization of expression data in a browser such as [IGV](http://software.broadinstitute.org/software/igv/). If no reference genome is available, one must supply RSEM a fasta file of transcript sequences. In addition, one can supply information that groups transcripts by gene, such that gene-level expression estimates.  

# Preliminaries

### Choosing an aligner
**STAR** is a splice-aware aligner that will produce gapped alignments to handle reads that span exon-intron junctions. Thus it is only appropriate when an annotated reference genome is available. STAR performs two-pass iterative mapping that identifies novel splice sites, and uses the updated splicing information to generate transcriptome (as well as genomic) alignments.

If one does not have a reference genome, **bowtie** and **bowtie2** are the other available aligner options. Note that, if one has an annotated genome and prefers mapping directly to the transcript sequences, one can also use these aligners. For RSEM to work properly when estimating expression from alignments to transcript sequences, those alignments must be ungapped. Bowtie is an ungapped aligner. While bowtie2 default behavior is to do gapped alignments, RSEM implements specific command line arguments such that it is run in an ungapped fashion. It has been shown that using bowtie2 this way is slightly more sensitive than bowtie, so we recommend it's use over bowtie2 unless there are project-specific reasons for using the former.

One has the ability to alter the command line arguments RSEM feeds to aligners, but one must be careful that those arguments don't produce alignments that RSEM cannot properly process, e.g. gapped alignments. 

### Gene vs. isoform level expression

RSEM has the ability to produce both gene and isoform-level expression estimates. However, accurate isoform level expression is typically much more challenging than gene-level estimation, and isoform-level estimates are far noisier. Thus, it is valuable to be able to group transcripts into genes. When the RSEM and aligner indices are built with an annotated reference genome, RSEM will automatically identify the genes to which annotated transcripts belong. Otherwise, one needs to provide a file that indicates the gene-transcript relationships. Depending upon the source of the transcriptome sequences and how they were annotated, this may or may not be straightforward (see below for an example).  


# Running RSEM

#### 1  Load the software and start SLURM session

RSEM is already installed as a module on the cluster. It can be loaded as follows:

	:::bash
	$ source new-modules.sh
	$ module load rsem/1.2.29-fasrc02

This loads RSEM and sequence alignment software, in this case Bowtie and Bowtie2.

A SAMtools module is also needed (SAMtools is software used to manipulate the SAM/BAM outputs of sequence alignment programs):

	:::bash
	$ module load samtools/1.4-fasrc01


For testing scripts and command line arguments with small data sets, it is possible to run RSEM in an interactive session. In all other cases, the compute resources and time required for running RSEM analyses will require running those analyses as SLURM jobs on the cluster.

When logged in on Odyssey, an interactive session can be launched with a command like the following, in which we are requesting one CPU (parameter "-n") in the interactive partition (parameter "-p"), for a bash session of 150 minutes ("-t") with memory 5000 Mb ("--mem="):

	:::bash
	$ srun -p interact -n 1 --pty -t 150 --mem=5000 /bin/bash


Create a separate working folder for this example, for data and results files, running jobs, etc., preferably on the /n/regal file system, which is recommended for SLURM jobs, especially if there is a lot of I/O or large files are written out during jobs.
For example, if you create your working directory under regal and name it "rsem_example", the full path to the directory could be:

	:::bash
	$ export RSEM_DIR=~/regal/rsem_example



#### 2  Download RNA-seq data


If you do not have you own RNA-seq fastq read data, you can download the following for this example.

The RNA-seq example data are from an experiment using human fibroblast cells, and are available on GEO at accession [GSE37704](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37704); there are six samples (three biological replicates in each of two conditions).
The six RNA-seq datasets (paired-end data in this example) can be downloaded from the GEO Short Read Archive as .sra datafiles; the SRA IDs of the datasets are "SRR493366", "SRR493367", "SRR493368" "SRR493369", "SRR493370", "SRR493371":

	:::bash
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145662/SRR493366/SRR493366.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145663/SRR493367/SRR493367.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145664/SRR493368/SRR493368.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145665/SRR493369/SRR493369.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145666/SRR493370/SRR493370.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145667/SRR493371/SRR493371.sra

These altogether can take about 15 mins to download to the cluster.

To extract fastq files from the GEO data, use NCBI program toolkit. It is already installed on Odyssey and can be loaded with:

	:::bash
	$ source new-modules.sh
	$ module load sratoolkit/0.2.3.4-2.bib-fasrc02

The fastq-dump command produces two compressed fastq files for each dataset, for the forward and reverse reads. In this example, this is done with a SLURM job-array:

	:::bash
	#!/bin/bash
	#SBATCH -n 1
	#SBATCH -N 1
	#SBATCH --mem 1000
	#SBATCH -p serial_requeue
	#SBATCH -o sra_%A_%a.out
	#SBATCH -e sra_%A_%a.err
	#SBATCH -J sra_arr
	#SBATCH -t 150
	#SBATCH --array=66-71

	source new-modules.sh
	module load sratoolkit/0.2.3.4-2.bib-fasrc02

	RSEM_DIR=~/regal/rsem_example

	srafile=SRR4933${SLURM_ARRAY_TASK_ID}.sra

	fastq-dump --split-files -I --gzip  -O $RSEM_DIR  $RSEM_DIR/$srafile


The fastq-dump processing of these datasets can take a couple of hours.  If you have much smaller input files, you may not have to use a SLURM job, but you could run on the interactive SLURM session.



#### 3  Transcriptome data


For transcriptome reference, we are using human transcriptome data from a recent annotation release of Ensembl. This dataset can be downloaded with:

	:::bash
	$ wget ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
	$ gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

The transcriptome data can be placed in a separate folder (the /n/regal file system is recommended), which we can mark with an environment variable, for easy reference:  

	:::bash
	$ export TRANS_DATA=/path/to/trasncriptome/data


A tab-delimited text file indicating the Gene/Transcript relationships (format: "gene_id(tab)transcript_id"), is also needed for RSEM. This file can be created using information from the transcriptome data source. In our example we obtained the table using sequence header annotation in the transcriptome fasta file.  
We have named this file "gene_transcript_map.txt". 

The human transcriptome data can be copied from the Informatics reference genome directory on Odyssey: /n/regal/informatics_public/ref/ensembl/release-83/homo_sapiens. This directory (/n/regal/informatics_public/ref) contains annotation datasets from other reference organisms that could potentially also be used as inputs for this RSEM example workflow. 




#### 4  Run RSEM

First process the reference transcriptome fasta file downloaded to build various index and other auxiliary files that RSEM needs; this is done with the RSEM "rsem-prepare-reference" function, which we run here as a SLURM job with the following script:

	:::bash
	#!/bin/bash
	#SBATCH -n 1
	#SBATCH -N 1
	#SBATCH --mem 8000
	#SBATCH -p serial_requeue
	#SBATCH -o rsem_prep.out
	#SBATCH -e rsem_prep.err
	#SBATCH -t 90

	source new-modules.sh
	module load rsem/1.2.29-fasrc02
	module load samtools/1.4-fasrc01

	TRANS_DATA=/path/to/trasncriptome/data

	rsem-prepare-reference --transcript-to-gene-map $TRANS_DATA/gene_transcript_map.txt --bowtie2 $TRANS_DATA/Homo_sapiens.GRCh38.cdna.all.fa  $TRANS_DATA/human38_cdna.rsem

Preparation of index files for this transcriptome took about 30 minutes.


Now we can quantify abundances of the transcripts in the RNA-seq datasets with the RSEM "rsem-calculate-expression" function. In the "rsem-calculate-expression" function we specify: 

- sequence alignment program to use, in this example Bowtie2 (parameter "--bowtie2"); 
- number of threads to use (parameter "--num-threads"); 
- paired-end fastq input files (parameter "--paired-end");
- and finally, as arguments, the paired-end fastq input datasets, the transcriptome index files, and the output prefix for each dataset. Note: the fastq inputs must be uncompressed for this example.

We run this as a SLURM job-array:

	:::bash
	#!/bin/bash
	#SBATCH -n 32
	#SBATCH -N 1
	#SBATCH --mem 75000
	#SBATCH -p general
	#SBATCH -o rsem_%A_%a.out
	#SBATCH -e rsem_%A_%a.err
	#SBATCH -J rsem_arr
	#SBATCH -t 1000
	#SBATCH --array=66-71


	source new-modules.sh
	module load rsem/1.2.16-fasrc03
	module load samtools/1.2-fasrc01

	TRANS_DATA=/path/to/trasncriptome/data
	RSEM_DIR=~/regal/rsem_example


	dataset=SRR4933${SLURM_ARRAY_TASK_ID}

	rsem-calculate-expression --bowtie2 --num-threads 32  --paired-end  $RSEM_DIR/${dataset}_1.fastq  $RSEM_DIR/${dataset}_2.fastq  $TRANS_DATA/human38_cdna.rsem  $RSEM_DIR/${dataset}.rsem


RSEM expression quantification processing of these datasets takes up to about 9-10 hours.



#### 5  RSEM output


The results of a RSEM run are written out in files that have the prefix indicated in the "rsem-calculate-expression" command. 

The main outputs of interest are the files containing the quantification results at the gene and isoform levels. The names of these files for each dataset are of the form:

- dataset_prefix.isoforms.results
- dataset_prefix.genes.results

These are tab-delimited files and contain expression estimates for each isoform ("transcript_id") or gene ("gene_id") as "expected_count", and also as TPM (Transcripts Per Million) and FPKM (Fragments Per Kilobase of transcript per Million mapped reads) numbers.



#### 6 References

Li, B. and Dewey, C. N. (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 12:323.
