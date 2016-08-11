Title: Kallisto example on Odyssey
Date: 2015-12-15 00:00
Category: Software
Tags: Next-Gen Sequencing, Transcriptome, RNA-seq Quantitation, Differential Expression, Kallisto
Summary: An example of quantifying RNA-seq expression with Kallisto on Odyssey cluster


[Kallisto](http://pachterlab.github.io/kallisto/) is a program for quantifying abundances of transcripts from RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads. It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment. The speed of kallisto makes it possible to use the bootstrap to determine uncertainty of estimates.


#### 1  Load the software and start SLURM session

Kallisto is already installed as a module on the cluster. It can be loaded as follows:

	:::bash
	$ source new-modules.sh
	$ module load gcc openmpi kallisto

To confirm that the module is loaded and check the version do:

	:::bash
	$ kallisto version

The result of the command should look like:

	kallisto, version 0.42.3


We recommend running longer kallisto jobs in an interactive SLURM session on Odyssey, or submit them as SLURM jobs.

When logged in on Odyssey, an interactive session can be launched with a command like the following, in which we are requesting one CPU (parameter "-n") in the interactive partition (parameter "-p"), for a bash session of 150 minutes ("-t") with memory 5000 Mb ("--mem="):

	:::bash
	$ srun -p interact -n 1 --pty -t 150 --mem=5000 /bin/bash


Create a separate working folder for this example, for data and results files, running jobs, etc., preferably on the /n/regal file system, which is recommended for SLURM jobs, especially if there is a lot of I/O or large files are written out during jobs.
For example, if you create your working directory under regal and name it "kallisto_example", the full path to the directory could be:

	:::bash
	$ export KALLISTO_DIR=~/regal/kallisto_example


#### 2  Download example data


If you don't have you own RNA-seq fastq read data and transcriptome reference fasta data, you can download the following for this example.

The RNA-seq example data are from an experiment using human fibroblast cells, and are available on GEO at accession [GSE37704](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37704); there are six samples (three biological replicates in each of two conditions).
The six RNA-seq datasets (paired-end data in this example) can be downloaded from the GEO Short Read Archive as .sra datafiles; the SRA ID's of the datasets are "SRR493366", "SRR493367", "SRR493368" "SRR493369", "SRR493370", "SRR493371":

	:::bash
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145662/SRR493366/SRR493366.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145663/SRR493367/SRR493367.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145664/SRR493368/SRR493368.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145665/SRR493369/SRR493369.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145666/SRR493370/SRR493370.sra
	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX145/SRX145667/SRR493371/SRR493371.sra

These altogether can take about 15 mins to download to the cluster.

To extract fastq files from the GEO data, use NCBI's program toolkit. It's already installed on Odyssey and can be loaded with:

	:::bash
	$ source new-modules.sh
	$ module load sratoolkit/0.2.3.4-2.bib-fasrc02

The fastq-dump command produces two compressed fastq files for each dataset, for the forward and reverese reads. In this example, this is done with a SLURM job-array:

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

	KALLISTO_DIR=~/regal/kallisto_example

	srafile=SRR4933${SLURM_ARRAY_TASK_ID}.sra

	fastq-dump --split-files -I --gzip  -O $KALLISTO_DIR  $KALLISTO_DIR/$srafile


The fastq-dump processing of these datasets can take a couple of hours.  If you have much smaller input files, you may not have to use a SLURM job, but you could run on the interactive SLURM session.


For human transcriptome reference, we are using transcriptome data from a recent annotation release of Ensembl. This dataset can be downloaded with:

	:::bash
	$ wget ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

The transcriptome data can be placed in a separate folder (the /n/regal file system is recommended), which we can mark with an environment variable, for easy reference:  

	:::bash
	$ export TRANS_DATA=/path/to/trasncriptome/data



#### 3  Run kallisto

First process the transcriptome fasta file downloaded to build a transcriptome index that kallisto needs with the kallisto "index" function; we name the index file "human_transcripts.idx":

	:::bash
	$ kallisto index -i $TRANS_DATA/human_transcripts.idx $TRANS_DATA/Homo_sapiens.GRCh38.cdna.all.fa.gz

Creation of the index takes about 10-15 minutes.


Now we can quantify abundances of the transcripts in the RNA-seq datasets with the kallisto "quant" function; kallisto can read in either plain-text or gzip-compressed files. In the "quant" function we specify: 

- location of the trasncriptome index (parameter "-i"); 
- kallisto output folder for each dataset (parameter "-o"); 
- number of bootstrap samples (parameter "-b"); 
- number of threads to use for bootstraping (parameter "-t"); 
- and finally, the location of the paired-end fastq, compressed or uncompressed, datasets (as arguments). 

We run this as a SLURM job-array:

	
	#!/bin/bash
	#SBATCH -n 32
	#SBATCH -N 1
	#SBATCH --mem 15000
	#SBATCH -p serial_requeue
	#SBATCH -o kallisto_%A_%a.out
	#SBATCH -e kallisto_%A_%a.err
	#SBATCH -J kallisto_arr
	#SBATCH -t 120
	#SBATCH --array=66-71


	source new-modules.sh
	module load gcc openmpi kallisto

	TRANS_DATA=/path/to/trasncriptome/data
	KALLISTO_DIR=~/regal/kallisto_example


	dataset=SRR4933${SLURM_ARRAY_TASK_ID}

	kallisto quant -i $TRANS_DATA/human_transcripts.idx -o $KALLISTO_DIR/${dataset}.kallisto_out -b 100 -t 32  $KALLISTO_DIR/${dataset}_1.fastq.gz  $KALLISTO_DIR/${dataset}_2.fastq.gz


Kallisto processing of these datasets takes about 30 minutes.


#### 4  Kallisto output


The results of a kallisto run are written in the specified output folder for each dataset (the -o option). The contents of each output folder are three files:

- abundance.tsv
- abundance.h5
- run_info.json

The main output of the kallisto quantification are the abundance estimates in the "abundance.tsv" file. For each transcript ("target_id"), abundance is reported in “estimated counts” (est_counts) and in Transcripts Per Million (tpm). 

The "run_info.json" file contains a summary of the run, including data on the number targets used for quantification, the number of bootstraps performed, the version of the program used and how it was called. 

The "abundance.h5" file contains the main quantification together with the boostraps in HDF5 compressed format, for large output of runs with many bootstraps.


#### 5 References

N. Bray, H. Pimentel, P. Melsted, and L. Pachter (2015). Near-optimal RNA-Seq quantification. [arXiv:1505.02710](http://arxiv.org/abs/1505.02710)

