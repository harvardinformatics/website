Title: RSEM example on Odyssey
Author: Adam Freedman
Date: 2017-06-19 00:00
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

RSEM has the ability to produce both gene and isoform-level expression estimates. However, accurate isoform level expression is typically much more challenging than gene-level estimation, and isoform-level estimates are far noisier. Thus, it is valuable to be able to group transcripts into genes. The grouping of transcripts into genes is carried out when the rsem (and associated aligner index) is built. If a genome sequence and an associated annotation file in either gtf or gff3 format are available for the species of interest, gene-transcript relationships can be automatically detected from the gtf file. In the absence of an annotated genome, a tab-separated text file can be supplied that defines those relationships: the first and second columns are gene id and transcript id, respectively, separated by a tab. Such an approach would be used for transcriptome assemblies where gene-isoform relationships are determined by the assembler, or where an annotation pipeline has been used to group putative transcript contigs into groups originating from the same gene. See below for examples of building the rsem index with each method.

### Pre-processing RNA-seq reads
Prior to conducting expression analyses, RNA-seq reads will have to be pre-processed. An important first step is to run [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to get a first glimpse at the quality of the data. Besides looking at base qualities, it is important to assess the extent to which there is residual adapter contamination after demultplexing. Depending upon the application and the type of sequencing library, one can use tools such as [TrimGalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to remove residual adapter and low quality bases. An example for running TrimGalore! on the odyssey cluster in the context of de novo transcriptome assembly can be found at [Trinity Best Practices](http://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html). It should be noted that the degree to which one should trim low quality bases will vary with the study. Current studies are not unanimous in the importance of quality trimming for de novo transcriptome assembly. For expression analyses, a recent study has suggested that quality-based trimming biases expression estimates. If you have questions regarding your particular study, please set up an appointment for a consultation with Informatics Group bioinformatics staff. 

# Running RSEM

#### 1  Loading RSEM and dependencies

RSEM is already installed as a module on the cluster. It can be loaded as follows:

	:::bash
	$ source new-modules.sh
	$ module load rsem/1.2.29-fasrc02

This loads RSEM and by default, the sequence aligners Bowtie and Bowtie2. If you plan on using STAR or another supported aligner, you will need to load it's module separately.

#### 2  Testing in interactive mode

For testing scripts and command line arguments with small data sets, it is possible to run RSEM in an interactive session. In all other cases, the compute resources and time required for running RSEM analyses will require running those analyses as SLURM jobs on the cluster.

When logged in on Odyssey, an interactive session can be launched with a command like the following, in which we are requesting one CPU (parameter "-n") in the interactive partition (parameter "-p"), for a bash session of 150 minutes ("-t") with memory 5000 Mb ("--mem="):

	:::bash
	$ srun -p interact -n 1 --pty -t 150 --mem=5000 /bin/bash


#Create a separate working folder for this example, for data and results files, running jobs, etc., preferably on the /n/regal file system, which is recommended for SLURM jobs, especially if there is a lot of I/O or large files are written out during jobs.
#For example, if you create your working directory under regal and name it "rsem_example", the full path to the directory could be:
#
#	:::bash
#	$ export RSEM_DIR=~/regal/rsem_example


#### 3a  Building an rsem index for an annotated reference genome
Assuming your genome fasta and annotation file are called genome.fasta and annotation.gtf, respectively, to build an rsem index for an annotated genome (that has a gtf or gff3 annotation file), using bowtie2 as the short read aligner, one would build an index as follows:

    :::bash
    #!/bin/bash
    #SBATCH -n 1
    #SBATCH --mem 8000
    #SBATCH -p serial_requeue
    #SBATCH -o rsem_prep.out
    #SBATCH -e rsem_prep.err
    #SBATCH -t 90

    source new-modules.sh
    module purge
    module load rsem/1.2.29-fasrc02

    TRANS_DATA=/path/to/trasncriptome/data

    rsem-prepare-reference --gtf $TRANS_DATA/annotation.gtf --bowtie2 $TRANS_DATA/genome.fasta  $TRANS_DATA/genome

This script tells rsem to 1) use annotation.gtf to determine isoform-gene relationships, build the aligner index files for bowtie2 for genome.fasta, and place the indices in the$TRANS_DATA directory with the prefix "genome". To use a gff3 annotation file instead of a gtf, you would use the flag --gff3 instead of --gtf. When either flag is observed, RSEM assumes that the fasta sequence being supplied is for a genome.

NOTE: UCSC gene annotation files downloaded as gtf do not preserve gene-isoform relationships. In this and other cases where gene and isoform ids are not linked in the annotation file, one can explicitly supply a tab-delimited text file indicating the Gene/Transcript relationships (format: "gene_id(tab)transcript_id"). In this scenario, one would build the RSEM index as follows:

    :::bash
    #!/bin/bash
    #SBATCH -n 1
    #SBATCH --mem 8000
    #SBATCH -p serial_requeue
    #SBATCH -o rsem_prep.out
    #SBATCH -e rsem_prep.err
    #SBATCH -t 90

    source new-modules.sh
    module purge
    module load rsem/1.2.29-fasrc02

    TRANS_DATA=/path/to/trasncriptome/data

    rsem-prepare-reference --gtf $TRANS_DATA/annotation.gtf --transcript-to-gene-map $TRANS_DATA/gene_transcript_map.txt --bowtie2 $TRANS_DATA/genome.fasta  $TRANS_DATA/genome  

NOTE: Even if you have a gtf file that properly links isoform and gene IDs, you can still override it by supplying a mapping file along with the --transcript-to-gene-map flag.  
	
#### 3b  Building an RSEM index for a set of transcript sequences
In cases where no genome and associated annotation file is available, one can build an index for the set of transcript sequences, supplying a transripts.fasta file instead of a genome fasta file. For Trinity de novo transcriptome assemblies, a transcript-to-gene-map can be easily created with RSEM's extract-transcript-to-gene-map-from-trinity tool. Alternatively, you can supply your own map file, e.g. based upon contig-gene relationships obtained from an annotation pipeline. In either case, one would build the index as follows:

    :::bash
    #!/bin/bash
    #SBATCH -n 1
    #SBATCH --mem 8000
    #SBATCH -p serial_requeue
    #SBATCH -o rsem_prep.out
    #SBATCH -e rsem_prep.err
    #SBATCH -t 03:00:00

    source new-modules.sh
    module purge
    module load rsem/1.2.29-fasrc02

    TRANS_DATA=/path/to/trasncriptome/data

    rsem-prepare-reference --transcript-to-gene-map $TRANS_DATA/gene_transcript_map.txt --bowtie2 $TRANS_DATA/transcripts.fasta  $TRANS_DATA/transcripts


#### 4  Estimate expression with RSEM
One then quantifies abundances of the transcripts in the RNA-seq dataset with the RSEM "rsem-calculate-expression" function. In the "rsem-calculate-expression" function we specify: 

- sequence alignment program to use, in this example Bowtie2 (parameter "--bowtie2"); 
- number of threads to use (parameter "--num-threads"); 
- paired-end fastq input files (parameter "--paired-end");
- and finally, as arguments, the paired-end fastq input datasets, the transcriptome index files, and the output prefix for each dataset. RSEM can handle either uncompressed or compressed (*gz) fastq files.

An example script for processing a single paired-end read sample, Chupacabra_R1.fq.gz and Chupacabra_R2.fq.gz is as follows

    :::bash
    #!/bin/bash
    #SBATCH -n 12
    #SBATCH -N 1
    #SBATCH --mem 12000
    #SBATCH -p serial_requeue
    #SBATCH -o rsem_%A.out
    #SBATCH -e rsem_%A.err
    #SBATCH -J rsem_chupa
    #SBATCH -t 10:00:00

    source new-modules.sh
    module load rsem/1.2.16-fasrc03
    module load samtools/1.2-fasrc01

    TRANS_DATA=/path/to/trasncriptome/data
    RSEM_DIR=~/regal/where/to/find/your/fastq/files

    rsem-calculate-expression --bowtie2 --num-threads 12 --paired-end  $RSEM_DIR/Chupacabra_R1.fq.gz  $RSEM_DIR/Chupacabra_R2.fq.gz  $TRANS_DATA/transcripts  $RSEM_DIR/Chupacabra_RSEM_bowtie2


It is easy to embed the above into either a bash script that loops over read pairs in order to quickly process the job submissions. Alternatively, it can be modified such that one supplies command line arguments for the fastq file and transcriptome locations.

Important things to note:
- For strand-specific RNA-seq protocols there are options to specify strandedness during the mapping phase. Consult your library kit manual, RSEM and short read aligner documentation to determine the approriate settings. In some cases, it is not immediately clear and requires some investigation!
- Depending upon the genome size/transcriptome complexity and the number of input reads, the amount of memory and time required to complete an RSEM job may vary. We recommend doing a test run with one data set, then using the SLURM sacct tool to get information on elapsed time and amount of memory used to optimize subsequent resource requests for subsequent job submissions.

#### 5  RSEM output

The results of a RSEM run are written out in files that have the prefix indicated in the "rsem-calculate-expression" command. 

The main outputs of interest are the files containing the quantification results at the gene and isoform levels. The names of these files for each dataset are of the form:

- dataset_prefix.isoforms.results
- dataset_prefix.genes.results

These are tab-delimited files and contain expression estimates for each isoform ("transcript_id") or gene ("gene_id") as "expected_count", and also as TPM (Transcripts Per Million) and FPKM (Fragments Per Kilobase of transcript per Million mapped reads) numbers.

#### 6 References

Li, B. and Dewey, C. N. (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 12:323.
