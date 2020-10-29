Title: RSEM example on Odyssey
Author: Adam Freedman
Date: 2020-10-28 00:00
Category: Tutorials
Tags: Next-Gen Sequencing, Transcriptome, RNA-seq Quantitation, Differential Expression, RSEM
Summary: An example of quantifying RNA-seq expression with RSEM on Odyssey cluster


[RSEM](http://deweylab.github.io/RSEM/README.html) is a software package for estimating gene and isoform expression levels from single-end or paired-end RNA-Seq data. The software works with transcriptome sequences and does not require a reference genome. It can either perform the read alignment step prior to quantification, or take an alignment (bam) file as input, so long as the alignment settings are appropriate for RSEM. Currently, RSEM can perform the alignment step with three different aligners: bowtie, bowtie2, or STAR. It uses the Expectation Maximization (EM) algorithm to estimate abundances at the isoform and gene levels.

If an annotated reference genome is available, RSEM can use a gtf file representation of those annotations to extract the transcript sequences for which quantification will be performed, and build the relevant genome and transcriptome indices. If this is done and the alignment step is implemented within RSEM, the option is available to also write the read alignments in genomic coordiinates, permitting visualization of expression data in a browser such as [IGV](http://software.broadinstitute.org/software/igv/). If no reference genome is available, one must supply RSEM a fasta file of transcript sequences. In addition, one can supply information that groups transcripts by gene, such that gene-level expression estimates. 

The example workflows below are demonstrated as sbatch scripts for using with the SLURM job scheduling engine. With modest tweaking,these can be reconfigured for other schedulers, e.g. LSF, SGE. 

### Choosing an aligner
[STAR](https://github.com/alexdobin/STAR)is a splice-aware aligner that will produce gapped alignments to handle reads that span exon-intron junctions. Thus it is only appropriate when an annotated reference genome is available. STAR performs two-pass iterative mapping that identifies novel splice sites, and uses the updated splicing information to generate transcriptome (as well as genomic) alignments.

If one does not have a reference genome, **bowtie** and **bowtie2** are the other available aligner options. Note that, if one has an annotated genome and prefers mapping directly to the transcript sequences, one can also use these aligners. For RSEM to work properly when estimating expression from alignments to transcript sequences, those alignments must be ungapped. Bowtie is an ungapped aligner. While bowtie2 default behavior is to do gapped alignments, RSEM implements specific command line arguments such that it is run in an ungapped fashion. It has been shown that using bowtie2 this way is slightly more sensitive than bowtie, so we recommend it's use over bowtie unless there are project-specific reasons for using the former.

One has the ability to alter the command line arguments RSEM feeds to aligners, but one must be careful that those arguments don't produce alignments that RSEM cannot properly process, e.g. gapped alignments. 

### Gene vs. isoform level expression

RSEM has the ability to produce both gene and isoform-level expression estimates. However, accurate isoform level expression is typically much more challenging than gene-level estimation, and isoform-level estimates are far noisier. Thus, it is valuable to be able to group transcripts into genes. The grouping of transcripts into genes is carried out when the rsem (and associated aligner index) is built. If a genome sequence and an associated annotation file in either gtf or gff3 format are available for the species of interest, gene-transcript relationships can be automatically detected from the gtf file. In the absence of an annotated genome, a tab-separated text file can be supplied that defines those relationships: the first and second columns are gene id and transcript id, respectively, separated by a tab. Such an approach would be used for transcriptome assemblies where gene-isoform relationships are determined by the assembler, or where an annotation pipeline has been used to group putative transcript contigs into groups originating from the same gene. See below for examples of building the rsem index with each method. Explicitly defining gene-transcript relationships may also be required if the gtf/gff3 annotation file does not have fields structured in a manner that RSEM is expecting. Note that a gtf/gff3 that does not include transcript entries but only includes gene-specific sub-features such as "exon" and "CDS" will not work for RSEM. The whole premise of RSEM is that modeling ambiguous read mappings across isoforms yields more robust gene-level expression estimates, such that it requires isoform/transcript annotation information.

### Pre-processing RNA-seq reads
Prior to conducting expression analyses, RNA-seq reads will have to be pre-processed. An important first step is to run [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to get a first glimpse at the quality of the data. Besides looking at base qualities, it is important to assess the extent to which there is residual adapter contamination after demultplexing. Depending upon the application and the type of sequencing library, one can use tools such as [TrimGalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to remove residual adapter and low quality bases. An example for running TrimGalore! on the odyssey cluster in the context of de novo transcriptome assembly can be found at [Trinity Best Practices](http://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html). It should be noted that the degree to which one should trim low quality bases will vary with the study. Current studies are not unanimous in the importance of quality trimming for de novo transcriptome assembly. For expression analyses, a recent study has suggested that quality-based trimming biases expression estimates. If you have questions regarding your particular study, please set up an appointment for a consultation with Informatics Group bioinformatics staff. 

### Workflows
While STAR can be run from within RSEM, this prevents best practice with respect to 2-pass mapping. In the 2-pass approach, one does a first round of alignment to collect sample-specific slice sites. In the second round, for each library one uses all of the splice site information obtained for all of the samples in the experiment. Thus, we present three different workflows:

* Using a reference genome, alignment with STAR outside of RSEM, with bam files provided to RSEM for abundance calculations (steps 1,2,3,and 6)
* Using a reference genome, wrapping of aligment step using bowtie2 from within RSEM (steps 4 and 7).
* Using a de novo transcriptome, wrapping of alignment step using bowtie2, while telling RSEM which contigs belong to which genes (steps 5 and 7).

We define these workflows as combinations of numbered steps enumerated below.

### Workflow steps (to be combined as described above)
Our example scripts used to implement workflows rely on creating an Anaconda environment, with doing so on the Cannon cluster described [here](https://informatics.fas.harvard.edu/python-on-cannon.html). Anaconda python distributions, available on the cluster, are used to create such an environment. To create an enviornment for rsem, and aligners it might use, one would do this:
    :::bash
    module load python
    conda create -n rsem -c bioconda rsem star bowtie2

You will likely get asked it you want to update particular dependencies, so click yes. Once you have created a conda environment, you can re-use it at any time. the -n specifies the name of the environment you care creating, the -c the channels you use to grab packages, followed by a list of packages you want installed into the environment. You can also add additional packages later to an environment. Examples for basic conda environment operations can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Scripts below demonstrate how to activate the enviroment.

#### 1. Build STAR index
        :::bash
        #!/bin/sh
        #SBATCH -n 12 
        #SBATCH -N 1
        #SBATCH -t 12:00:00
        #SBATCH --mem=64000
        #SBATCH -p shared
        #SBATCH -e starprep.e
        #SBATCH -o starprep.o

        module purge
        module load python
        source activate rsem

        # $1 == 1 - read length, i.e. if you did 2 x 100 PE, this value is 99
        # $2 == gtf annotation file
        # $3 == full path to genome fasta
        # Note: this script writes the STAR index to the current working directory where
        # the SLURM script is being executed.

        STAR --runMode genomeGenerate -runThreadN 12 --sjdbOverhang $1 --sjdbGTFfile $2 --genomeDir $(pwd) --genomeFastaFiles $3 

Note: building a STAR index can be a memory-itensive process, and one may need to allocate more memory to the job. If the SLURM jobs fails to exceeding the memory allocation, the error log will indicate the minimum amount of memory one should reserve.
#### 2. Do first-pass alignment to the genome with STAR
        :::bash
        #!/bin/bash 
        #SBATCH -N 1
        #SBATCH -n 16
        #SBATCH -p shared
        #SBATCH -e star_%A.e
        #SBATCH -o star_%A.o
        #SBATCH -J star
        #SBATCH --mem=64000  
        #SBATCH -t 23:00:00 

        module purge
        module load python
        source activate rsem
	
        # $1 == path to directory where STAR index lives
        # $2 == a prefix for the STAR output, typically includes sample name
        # $3 == R1 fastq file
        # $4 == R2 fastq file
        # The readFilesCommand is set to zcat under the assumption that the fastq files are provided
        # as gzipped files. The number of threads, total memory, and length of job will depend upon
        # a number of factors including the library size, size of genome, etc. 	

        STAR --runThreadN 16 --genomeDir $1 --readFilesCommand zcat --outFileNamePrefix $2 --readFilesIn $3 $4

#### 3. 2nd pass STAR alignments
As with 1st pass, these are done for each sample. The two differences are that 1) we use splice sites detected from all 1st-pass alignments performed as part of your experiment and 2) we tell STAR to output a bam in transcriptome coordinates (with *Aligned.toTranscriptome.out.bam* as a suffix), which is what one will then feed into RSEM.
        :::bash
        #!/bin/bash
        #SBATCH -N 1
        #SBATCH -n 16
        #SBATCH -p shared
        #SBATCH -e star_%A.e
        #SBATCH -o star_%A.o
        #SBATCH -J star
        #SBATCH --mem=64000
        #SBATCH -t 23:00:00

        module purge
        module load python
        source activate rsem

        # $1 == path to directory where STAR index lives
        # $2 == space separated list of all splice site *tab files generated from 1-st pass
        # $3 == a prefix for the STAR output, typically includes sample name
        # $4 == R1 fastq file
        # $5 == R2 fastq file
	
        STAR --runThreadN 16 --genomeDir $1 --sjdbFileChrStartEnd $2 --quantMode TranscriptomeSAM --readFilesCommand zcat --outFileNamePrefix $3 --readFilesIn $4 $5

#### 4. Build an RSEM index for wrapping bowtie2 alignment to genome annotations
        :::bash
        #!/bin/bash 
        #SBATCH -n 6
        #SBATCH –N 1
        #SBATCH -p shared
        #SBATCH -e rsemindex_%A.err
        #SBATCH -o rsemindex_%A.out
        #SBATCH –J rsemindex
        #SBATCH --mem=8000
        #SBATCH -t 03:00:00

        module purge
        module load python
        source activate rsem
	
        # $1 == path to gtf annotation file
        # $2 == genome fasta file
        # $3 == name of index, e.g. the genome name without the .fa/.fasta file extension
        # The --bowtie2 argument tells RSEM to build a bowtie2 index for the alignment step.
        # On Odyssey, loading the RSEM module also loads a compatible version of bowtie 2.

        rsem-prepare-reference -p 6 --bowtie2 --gtf $1 $(pwd)/${2} $(pwd)/$3 

#### 5. Build an RSEM index, specifying gene-transcript relationsips
This approach is used when gene-transcript relationships are defined from an external source,e.g. if you are using BLASTP/X annotation of de novo transcriptome assembly contigs to associated contigs to gene-level symbols.
        :::bash
        #!/bin/bash 
        #SBATCH -n 6
        #SBATCH –N 1
        #SBATCH -p shared
        #SBATCH -e rsemindex_%A.err
        #SBATCH -o rsemindex_%A.out
        #SBATCH –J rsemindex
        #SBATCH --mem=8000
        #SBATCH -t 03:00:00

        module purge
        module load python
        source activate rsem

        # $1 == path to gtf annotation file
        # $2 == tab-separated file with gene id and transcript/contig id in 1st and 2nd columns, respectively
        # $3 == genome fasta file
        # $4 == name of index, e.g. the genome name without the .fa/.fasta file extension
        # The --bowtie2 argument tells RSEM to build a bowtie2 index for the alignment step.
        # On Odyssey, loading the RSEM module also loads a compatible version of bowtie 2.

        rsem-prepare-reference -p 6 --bowtie2 --gtf $1 --transcript-to-gene-map $2 $(pwd)/${3} $(pwd)/$4

NOTE: UCSC gene annotation files downloaded as gtf do not preserve gene-isoform relationships. In this and other cases where gene and isoform ids are not linked in the annotation file, the above approach can be used. Even if you have a gtf file that properly links isoform and gene IDs, you can still override it by supplying a mapping file along with the --transcript-to-gene-map flag.  
	

### Estimating expression
One quantifies abundances of the transcripts in the RNA-seq dataset with the RSEM "rsem-calculate-expression" function. In the "rsem-calculate-expression" function we specify:

- sequence alignment program to use (parameter "--bowtie2", or "--star");
- number of threads to use (parameter "--num-threads");
- either arguments specifying fastq read files or a pre-computed alignment (bam) file, as described below;
- the name (including full path) to the RSEM index, and;
- the sample name

#### 6. Estimate expression from pre-computed STAR alignments  
        :::bash
        #!/bin/bash
        #SBATCH -n 16
        #SBATCH -N 1
        #SBATCH --mem 64000
        #SBATCH -p serial_requeue,shared
        #SBATCH -o rsem_%A.out
        #SBATCH -e rsem_%A.err
        #SBATCH -J rsem
        #SBATCH -t 10:00:00

        module purge
        module load python
        source activate rsem
	
        # $1 == STAR alignments bam file
        # $2 == name of STAR index
        # #3 == sample name

        rsem-calculate-expression --star --num-threads 16 --alignments $1 $(pwd)/$2 $(pwd)/$3 

#### 7. Align reads with bowtie2 and estimate expression
        :::bash
        #!/bin/bash
        #SBATCH -n 16
        #SBATCH -N 1
        #SBATCH --mem 64000
        #SBATCH -p serial_requeue,shared
        #SBATCH -o rsem_%A.out
        #SBATCH -e rsem_%A.err
        #SBATCH -J rsem
        #SBATCH -t 10:00:00

        module purge
        module load python
        source activate rsem

        # $1 == R1
        # $2 == R2
        # $3 == path to rsem index
        $ $4 sample name        
 
        rsem-calculate-expression --bowtie2 --num-threads 16 --paired-end $1 $2 $3 $4 

#### Important things to note:
- It is easy to embed the above into either a bash script that loops over read pairs in order to quickly process the job submissions. Alternatively, it can be modified such that one supplies command line arguments for the fastq file and transcriptome locations.

- For strand-specific RNA-seq protocols there are options to specify strandedness during the mapping phase. Consult your library kit manual, RSEM and short read aligner documentation to determine the approriate settings. Typically, for dUTP libraries, one uses "--forward-prob 0", and for ligation-stranded protocols, one uses "--forward-prob 1". In some cases, it is not immediately clear and requires some investigation!
- Depending upon the genome size/transcriptome complexity and the number of input reads, the amount of memory and time required to complete an RSEM job may vary. We recommend doing a test run with one data set, then using the SLURM sacct tool to get information on elapsed time and amount of memory used to optimize subsequent resource requests for subsequent job submissions.

#### Understanding RSEM output

The results of a RSEM run are written out in files that have the prefix indicated in the "rsem-calculate-expression" command. 

The main outputs of interest are the files containing the quantification results at the gene and isoform levels. The names of these files for each dataset are of the form:

- dataset_prefix.isoforms.results
- dataset_prefix.genes.results

These are tab-delimited files and contain expression estimates for each isoform ("transcript_id") or gene ("gene_id") as "expected_count", and also as TPM (Transcripts Per Million) and FPKM (Fragments Per Kilobase of transcript per Million mapped reads) numbers.

#### 6 References

Li, B. and Dewey, C. N. (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 12:323.
