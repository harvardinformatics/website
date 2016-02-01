Title: Transcriptome assembly and quantification with Trinity, Bowtie and RSEM
Date: 2013-12-07 00:00
Category: Software
Tags: RNA-Seq Analysis,De Novo Assembly,Trinity,Bowtie,RSEM
Summary: Assembly and quantification of transcriptome data with Trinity, Bowtie, and RSEM

The pipeline below can be used for _de novo_ assembly of transcripts using RNA-Seq read data. This assembly pipeline could be used, for instance, when there is no reference sequence to map reads to. First, reads from all RNA-Seq samples are merged, and this merged read set is used to assemble a transcriptome with Trinity. Next, reads of each sample are aligned to the transcriptome using Bowtie. Finally, abundance values are calculated with RSEM.To begin, the following Trinity and Bowtie versions were loaded:

    module load bio/trinityrnaseq_r2012-10-05
    module load bio/bowtie-0.12.7

The path to this Trinity installation on the cluster is at the following location:

    TRINITYDIR="/n/sw/trinityrnaseq_r2012-10-05"

Note that before using this pipeline, the read data may need to be trimmed and/or filtered. 


1. Run Trinity to generate a transcriptome assembly Here we run Trinity with a merged set of reads, which is especially appropriate if individual sample datasets are small.

        $TRINITYDIR/Trinity.pl --seqType fq --left all_samples.R1.fastq.gz --right all_samples.R2.fastq.gz --no_cleanup --JM 100G --CPU 8 --inchworm_cpu 8 --bflyHeapSpaceMax 20G --bflyCPU 4
    
    This creates, among other outputs, a FASTA file with the assembly sequences.  This file can optionally be filtered (e.g. by length) before further analysis:
    
        trinity_out_dir/Trinity.fasta
        
2. Align reads to the transcripts using Bowtie Run the following for each sample of interest (or grouped subsets of samples) against the assembly created in the first step (`Trinity.fasta`):
    
        $TRINITYDIR/util/alignReads.pl --seqType fq --left sampleX.R1.fastq --right sampleX.R2.fastq --target trinity_out_dir/Trinity.fasta --aligner bowtie -o sampleX.bowtie_output -- -p 16
    
    Bowtie output (BAM files, principally) will be located in folder below:
    
        sampleX.bowtie_output
        
3. Read quantitation with RSEM Run RSEM (found in Trinity installation) to quantify reads in transcripts. Use the assembly from step 1 (`Trinity.fasta`) and Bowtie output from step 2 (`sampleX.bowtie_output.nameSorted.bam`). Run the following command separately for each sample:
    
        $TRINITYDIR/util/RSEM_util/run_RSEM.pl --transcripts trinity_out_dir/Trinity.fasta --name_sorted_bam sampleX.bowtie_output/sampleX.bowtie_output.nameSorted.bam --paired
    
    RSEM generates tables for each sample containing read abundance values (including FPKM values), by gene and by transcript isoform:
    
        RSEM.genes.results
        RSEM.isoforms.results