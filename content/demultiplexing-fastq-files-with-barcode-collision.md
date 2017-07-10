Title: Demultiplexing FASTQ files with barcode collision
Date: 2014-06-02 00:00
Category: Tutorials
Tags: Next-Gen Sequencing
Summary: Multiplexing allows several samples to be sequenced in the same lane, but there can be problems with very short barcodes.  

Multiplexing allows several samples to be sequenced in the same lane. During library prep, each read is tagged with a barcode which serves as the sample label. Usually, the [hamming distance](http://en.wikipedia.org/wiki/Hamming_distance) between barcodes is high enough that, even with 1 or 2 mismatches allowed, there is no ambiguity about which sample a given read belongs to. Unfortunately, with very short barcodes, this is not always the case. Consider the following barcodes: 

    Sample_1   ACAGTG 
    Sample_2   CTTGTA 
    Sample_3   CGATGT 

With 2 mismatches, an index read of, say, CTAGTG maps equally well to Sample_1 and Sample_2. It is usually a good idea to know when this happens, so many demultiplexing tools will error when a read maps to more than one barcode. However, instead of reducing the number of mismatches, we may prefer to retain only those reads that map to exactly one barcode. The script [demultiplex.pl](https://github.com/harvardinformatics/bioblog/tree/master/Illumina.Demultiplex) performs this task. For example, to demultiplex R1.fq with 2 mismatches, use the following command:

    :::shell-session
    demultiplex.pl -fastq R1.fq -mismatches 2 -barcodes mybarcodes.txt


The barcodes file must be a tab-delimited, two-column file, where the first column contains sample names and the second column the corresponding barcodes. For a complete list of input options, type

    :::shell-session
    demultiplex.pl -h
