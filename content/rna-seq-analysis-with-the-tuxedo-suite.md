Title: RNA-Seq analysis with the Tuxedo suite 
Date: 2013-12-08 00:00
Category: Blog
Tags: RNA-Seq Analysis,cuffdiff,cufflinks,tophat,cuffmerge
Summary: An example of a pipeline for analyzing RNA-seq data with the Tuxedo suite, primarily for assessing differential expression between samples.

Here is an example of a pipeline for analyzing RNA-seq data with the Tuxedo suite, primarily for assessing differential expression between samples. In this example, we use two samples: sampleX and sampleY. The input reads (single-end) are in fastq files. The reference genome data is in files "genome.fa" (sequence) and "genes.gtf" (gene annotation). First, align reads with splice junction mapper TopHat, using Bowtie:

    :::shell-session
    $ tophat --bowtie1 -G genes.gtf --no-novel-juncs --transcriptome-only -o sampleX  reference/genome  sampleX.R1.fastq.gz
    $ tophat --bowtie1 -G genes.gtf --no-novel-juncs --transcriptome-only -o sampleY  reference/genome  sampleY.R1.fastq.gz

The aligned reads are output in "accepted_hits.bam" file. Next, assemble transcripts and estimate read abundances with Cufflinks:

    :::shell-session    
    $ cufflinks -p 8 -G genes.gtf  --multi-read-correct --frag-bias-correct genome.fa -o sampleX/cufflinks_out  sampleX/accepted_hits.bam
    $ cufflinks -p 8 -G genes.gtf  --multi-read-correct --frag-bias-correct genome.fa -o sampleY/cufflinks_out  sampleY/accepted_hits.bam

Merge assemblies with Cuffmerge:

    :::shell-session
    $ cuffmerge -p 8 -o cuffmerge_out -g genes.gtf -s genome.fa assembly_list.txt

The "assembly_list.txt" file contains a list of (full) paths of the "transcripts.gtf" output of Cufflinks, e.g.: 

    sampleX/cufflinks_out/transcripts.gtf 
    sampleY/cufflinks_out/transcripts.gtf 
    
Finally, find differences in transcript expression with Cuffdiff:

    :::shell-session
    $ cuffdiff -p 16 -o cuffdiff_out  -N -c 15 -b genome.fa -u -L sampleX,sampleY  cuffmerge_out/merged.gtf  sampleX/accepted_hits.bam  sampleY/accepted_hits.bam