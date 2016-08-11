Title: Analyzing public ChIP-seq data for orthologs of genes of interest
Date: 2016-03-10 00:00
Category: Tutorials
Tags: Next-Gen Sequencing, ChIP-Seq, Peaks, Orthologs 
Summary: An example workflow for identifying orthologs of genes of interest from one organism by analyzing ChIP-seq data from another organism


#### 1  Introduction

With all the next-generation sequencing data now becoming available from various publications and in public databases, like the [GEO](http://www.ncbi.nlm.nih.gov/geo) database at NCBI, it has been become increasingly useful to mine these data for re-analysis or meta-analysis in other research studies.
In this example, given a set of genes of interest in one organism, called here the "original" organism, for instance a group of co-regulated genes with established roles in a developmental process of that organism, we show how to analyze ChIP-seq data to search for orthologs in another organism, called here the "target" organism.

Public ChIP-seq data can be found in various formats, which would determine different starting points for data processing, for example: 

- A. Data can be in the form of read sequences, in FASTQ format from the original sequencing experiment output. In this case, reads would need to be mapped to the organism reference genome with some alignment program, for example [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and then peaks can be called from the alignmengt results with software like [MACS2](https://github.com/taoliu/MACS/); see section 3 below.

- B. Public data can also be in the form of alignment results, e.g. SAM/BAM files, or perhaps bigWig, bedGraph or similar format, in which case, they can be processed with peak-calling software like MACS2, after possibly some necessary reformatting; see section 4 below.

- C. Called ChIP-seq peaks may also be available for more than one replicates in the datasets of interest. In this case the replicate sets of peaks may need to be combined into one set, by taking intersections or unions, or combining the sets using additional information to filter peaks, e.g. peak confidence scores; see section 5 below.

- D. Ideally, the public data provide a single set of peaks for each condition/dataset of interest, with their positions listed in a table, like in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format, with chromosome, start and end information. This would be the end point of all the other pre-proprocessing in the cases described above. 



#### 2  Load software

This workflow assumes a Unix (or Linux) commandline setting. Other than standard Unix commands, some additional software tools are needed. 
Many software programs needed for ChIP-seq analysis are already installed on the Harvard FAS Odyssey cluster as modules and can be loaded using the Lmod system. To activate Lmod use:

	$ source new-modules.sh

Programs can then be loaded with the "module load" command. Here is how we can load a few programs needed below. 

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), a program for checking the quality of FASTQ read data, is installed as a module on the cluster. It can be loaded as follows:

	$ module load fastqc
	
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) trimming tool for cleaning and filtering Illumina NGS single-end and paired-end reads, can be easily installed locally, by downloading and unzipping the binary executable:

	$ wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip
	$ unzip Trimmomatic-0.35.zip

You can mark the location of the Trimmomatic directory with an environment variable, for easy reference:  
	
	$ export TRIMMOMATIC_HOME=/path/to/trimmomatic/Trimmomatic-0.35

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is a fast and memory-efficient software program for aligning sequencing reads to long reference sequences: 

	$ module load bowtie2

Further programs that are needed for this workflow will be introduced in sections below.

We recommend running longer jobs as SLURM jobs on the cluster, or in an interactive SLURM session on Odyssey, for smaller datasets.

When logged in on Odyssey, an interactive session can be launched with a command like the following, in which we are requesting one CPU (parameter "-n") in the interactive partition (parameter "-p"), for a bash session of 150 minutes ("-t") with memory 5000 Mb ("--mem="):
	
	$ srun -p interact -n 1 --pty -t 150 --mem=5000 /bin/bash

Create a separate working folder for this example, for data and results files, running jobs, etc., preferably on the /n/regal file system, which is recommended for SLURM jobs, especially if there is a lot of I/O or large files are written out during jobs.
For example, if you create your working directory under regal and name it "chip_analysis", the full path to the directory could be:

	$ export CHIP_DIR=~/regal/chip_analysis



#### 3 Analysis starting from read sequence data

If starting the ChIP-seq analysis from NGS sequencing data in the form of FASTQ files (see bullet A in Introduction), we will need to assess the quality of the reads, trim low quality reads (or remove any primers from reads), and map reads to the genome.

Let us suppose we start with single-end FASTQ data from a ChIP-seq replicate, with a treatment and a control file:

	rep1_treatment.fq.gz
	rep1_control.fq.gz


Check quality statistics of reads using FastQC. This is done here using a SLURM job:

	#!/bin/bash   
	#
	#SBATCH -n 2  
	#SBATCH -N 1  
	#SBATCH -t 50  
	#SBATCH -p serial_requeue 
	#SBATCH --mem=3000  
	#SBATCH -o fastqc.out  
	#SBATCH -e fastqc.err 
	#SBATCH -J fastqc

	source new-modules.sh
	module load fastqc

	CHIP_DIR=~/regal/chip_analysis

	mkdir $CHIP_DIR/fastqc_out

	fastqc -t 2 --casava --noextract --nogroup --outdir $CHIP_DIR/fastqc_out  \
	   $CHIP_DIR/rep1_treatment.fq.gz  $CHIP_DIR/rep1_control.fq.gz 

The FastQC results can be viewed in .html files which are found in the output folder we created: $CHIP_DIR/fastqc_out


If the quality stats show that there are many reads with low quality bases, or there are sequencing primers in the reads (or other such contaminating sequences), trim with Trimmomatic using the following SLURM job:

	#!/bin/bash
	#
	#SBATCH -n 8
	#SBATCH -N 1
	#SBATCH --mem=3000
	#SBATCH -p serial_requeue
	#SBATCH -o trim.out
	#SBATCH -e trim.err
	#SBATCH -J trim_job
	#SBATCH -t 50

	TRIMMOMATIC_HOME=/path/to/trimmomatic/Trimmomatic-0.35
	CHIP_DIR=~/regal/chip_analysis

	java -jar $TRIMMOMATIC_HOME/trimmomatic-0.35.jar SE -threads 8 -phred33  \
	 $CHIP_DIR/rep1_treatment.fq.gz $CHIP_DIR/rep1_treatment.trimmed.fq.gz  \
	 ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

	java -jar $TRIMMOMATIC_HOME/trimmomatic-0.35.jar SE -threads 8 -phred33  \
	 $CHIP_DIR/rep1_control.fq.gz $CHIP_DIR/rep1_control.trimmed.fq.gz  \
	 ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20


The cleaned reads in output files "rep1_treatment.trimmed.fq.gz" and "rep1_control.trimmed.fq.gz" can now be mapped with Bowtie2 to the reference genome, in this example the human reference genome version hg38.

This FAS Informatics directory "/n/regal/informatics_public/ref" contains reference genome sequences and indices as well as annotation datasets that can be used with Bowtie2 (and other programs).

The Bowtie2 index "/n/regal/informatics_public/ref/ucsc/Homo_sapiens/hg38/Sequence/Bowtie2Index" is used in the example script below:

If needed, other preprocessed reference data (several organisms and versions) can be downloaded from the Illumina [iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html) website.

	#!/bin/bash
	#
	#SBATCH -n 16
	#SBATCH -N 1
	#SBATCH --mem 15000
	#SBATCH -p serial_requeue
	#SBATCH -o bowtie2.out
	#SBATCH -e bowtie2.err
	#SBATCH -J bowtie2
	#SBATCH -t 1200
	
	source new-modules.sh
	module load bowtie2
	
	GENOME=/n/regal/informatics_public/ref/ucsc/Homo_sapiens/hg38/Sequence/Bowtie2Index
	CHIP_DIR=~/regal/chip_analysis

	bowtie2 -p 16 -x $GENOME/genome -U $CHIP_DIR/rep1_treatment.trimmed.fq.gz -S $CHIP_DIR/rep1_treatment.sam
	
	bowtie2 -p 16 -x $GENOME/genome -U $CHIP_DIR/rep1_control.trimmed.fq.gz -S $CHIP_DIR/rep1_control.sam


The aligned reads are in SAM output files "rep1_treatment.sam" and "rep1_control.sam".



#### 4 Analysis starting from mapped read sequences 

If starting with aligned reads in a format like SAM (see bullet B in Introduction), it would be good to sort the reads, if they are not already sorted. In this example, we are sorting reads by genome position coordinate using [Picard](http://broadinstitute.github.io/picard/) tools, which are installed on the cluster and can be loaded with "module load picard-tools".

	#!/bin/bash
	#
	#SBATCH -n 1
	#SBATCH -N 1
	#SBATCH --mem 20000
	#SBATCH -p serial_requeue
	#SBATCH -o sort.out
	#SBATCH -e sort.err
	#SBATCH -J sort
	#SBATCH -t 200
	
	source new-modules.sh
	module load picard-tools 
	
	CHIP_DIR=~/regal/chip_analysis
	
	java -Xmx5g -jar $PICARD_TOOLS_HOME/SortSam.jar  I=$CHIP_DIR/rep1_treatment.sam  \
	   O=$CHIP_DIR/rep1_treatment.sorted.bam  SO=coordinate

	java -Xmx5g -jar $PICARD_TOOLS_HOME/SortSam.jar  I=$CHIP_DIR/rep1_control.sam  \
	   O=$CHIP_DIR/rep1_control.sorted.bam  SO=coordinate

Having obtained the sorted read alignment files (in BAM format) "rep1_treatment.sorted.bam" and "rep1_control.sorted.bam", we can now proceed to call peaks from the data. 
We are using [MACS2](https://github.com/taoliu/MACS) to call peaks; MACS2 assesses genome complexity to evaluate the significance of enriched ChIP regions and identifies transcript factor binding sites. It is installed on the cluster and can be loaded with "module load macs2".
The macs2 command in the following SLURM script takes the BAM files as input and produces outputs prefixed with "rep1"; the effective genome size parameter is set for a human size reference genome ("-g hs").

	#!/bin/bash
	#
	#SBATCH -n 1
	#SBATCH -N 1
	#SBATCH --mem 15000
	#SBATCH -p serial_requeue
	#SBATCH -o macs2.out
	#SBATCH -e macs2.err
	#SBATCH -J macs2
	#SBATCH -t 600
	
	source new-modules.sh
	module load macs2
	
	CHIP_DIR=~/regal/chip_analysis
	
	macs2 callpeak -t $CHIP_DIR/rep1_treatment.sorted.bam -c $CHIP_DIR/rep1_control.sorted.bam -f BAM -B -g hs -n rep1

The main MACS2 output file that we are interested in for this workflow is the one with suffix "_peaks.narrowPeak"; it is a table listing information about the peaks, including location and significance value.

We can select and sort a subset of the columns of this table, "chr", "start", "end", "peak_ID", "fold_enrichment", "-log10(qvalue)", with:

	$ awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$7,$9}' rep1_peaks.narrowPeak | sort -k1,1d -k2,3n | uniq > rep1_peaks.bed

So now we have a sorted list of peaks in BED format "rep1_peaks.bed".



#### 5 Combining ChIP-seq peak results from replicates 

If ChIP-seq peaks are available from more than one replicates (see bullet C in Introduction), the replicate sets of peaks can be compared to see how well they agree and may be combined into one set. For instance, if we have peaks from two replicates in BED files, "rep1_peaks.bed" and "rep2_peaks.bed", we can look at their intersection; according to [ENCODE project guidelines](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/), more than two ChIP-seq replicates do not seem to significantly improve site discovery.

The main package needed to process data in BED and related formats is [bedtools](http://bedtools.readthedocs.org) (see [bedtools installation information](http://bedtools.readthedocs.org/en/latest/content/installation.html) for installation on a desktop).

If using the Odyssey cluster, the bedtools package is already installed and can be loaded as follows:

	$ module load bedtools2

We can take the intersection of the two sets of peaks (and sort the outcome) with:

	$ bedtools intersect -a rep1_peaks.bed -b rep2_peaks.bed | sort -k1,1d -k2,3n | uniq > target_chip_peaks.bed

The resulting file with the common peak intervals is "target_chip_peaks.bed".

As a simple rule of thumb from the ENCODE project, if 75% of peaks overlap between replicates, then this indicates good agreement between the replicates.


##### IDR

Peak replicate agreement can also be examined in a more statistically rigorous manner, by taking peak score rankings into account.
If two replicates measure the same underlying biology, the most significant peaks, which are likely to be genuine signals, are expected to have high consistency between replicates, whereas peaks with low significance, which are more likely to be noise, are expected to have low consistency. This consistency can be used to evaluate the significance of peaks, through the Irreproducible Discovery Rate [IDR](http://projecteuclid.org/download/pdfview_1/euclid.aoas/1318514284).


[IDR software](https://sites.google.com/site/anshulkundaje/projects/idr#TOC-Code-for-IDR-Analysis) can be downloaded and uncompressed with:

	$ wget https://sites.google.com/site/anshulkundaje/projects/idr/idrCode.tar.gz
	$ tar xvfz idrCode.tar.gz

Change directory to where the "idrCode" folder has been uncompressed, for instance: 

	$ cd /path/to/idrCode

and load R:

	$ module load R
 
IDR analysis can be run using, as input, MACS2 "_peaks.narrowPeak" output files with the following command, which specifies the prefix "common" for output files and sets the peak "p.value" as the ranking measure for the IDR analysis:

	$ Rscript batch-consistency-analysis.r $CHIP_DIR/rep1_peaks.narrowPeak $CHIP_DIR/rep2_peaks.narrowPeak -1 common 0 F p.value

One of the output files, "common-overlapped-peaks.txt", lists overlapping pairs of peaks from the two replicates with the their IDR level; here is an example of the first few lines of that file:

	"chr1" "start1" "stop1" "sig.value1" "chr2" "start2" "stop2" "sig.value2" "idr.local" "IDR"
	"1" "chr22" 211 504 13.3310166666667 "chr22" 179 497 17.1514933333333 0.21553976929997 0.0349000090705554
	"2" "chr3" 961 1189 29.07822 "chr3" 972 1194 24.23245 0.0209891604384735 0.00215639386575345
	"3" "chr3" 1254 1450 7.6089 "chr3" 1272 1429 7.90452 0.272450030080924 0.0495086468630693
	"4" "chr11" 2162 2430 14.350225 "chr11" 2162 2395 13.5125 0.78864899085754 0.148831691925851


Another output file, "common-npeaks-aboveIDR.txt", contains the number of overlapping pairs of peaks at different IDR levels. Here is an example of the contents of that file:

	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.01 33047 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.02 36755 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.03 39387 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.04 41641 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.05 43668 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.06 45494 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.07 47089 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.08 48468 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.09 49667 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.1 50786 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.11 51849 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.12 52865 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.13 53819 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.14 54725 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.15 55599 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.16 56468 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.17 57334 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.18 58193 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.19 59050 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.2 59911 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.21 60771 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.22 61639 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.23 62517 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.24 63406 
	rep1_peaks.narrowPeak rep2_peaks.narrowPeak 0.25 64307 


To get peaks below below a certain IDR level, e.g IDR less or equal to 0.10, filter entries in the "common-overlapped-peaks.txt" output file using the following command:

	$ awk '$11 <= 0.10 && NR > 1 {print $2"\t"$3"\t"$4"\t"$11}' common-overlapped-peaks.txt | sed -e 's/"//g' | sort -k1,1d -k2,3n | uniq > idr_0.10_peaks.bed


The resulting sorted peak list will be in "idr_0.10_peaks.bed".



#### 6  Ortholog analysis starting with a set of ChIP-seq peaks: gather datafiles needed

If the starting point is one set of ChIP-seq peaks (see bullet D in Introduction), we'll work on this peak list (in BED format) together with appropriate annotation data.

So the data files we'll need for searching for orthologs of interest based on ChIP-seq data are:

- list of ChIP-seq peaks in target organism, in BED format, as described above, named "target_chip_peaks.bed" in this workflow.
Example of peak list tab-delimited format (chromosome, start, end):

```
chr1	829619	830620
chr1	832973	833974
chr1	835512	836513
```

If the peak BED file has a few additional columns, for example a column for peak ID and/or peak score, it can still be processed by the workflow below.

- list of genes of interest in original organism, say worm, named "genes_of_interest.txt"
This table may optionally have additional columns, but only the first one is needed for the workflow below. Example of the data format (original gene name):

```
ajm-1
cdh-4
crb-1
src-1
```

Other annotation data needed:

- list of genes and their positions in target organism (e.g. human), named "target_gene_annotation.bed" (chromosome names sorted alphabetically)
Example of the tab-delimited data format (chromosome, start, end, gene ID):

```
chr1	11869	14409	ENSG00000223972
chr1	14404	29570	ENSG00000227232
chr1	17369	17436	ENSG00000278267
```

- list of target and original organism pairs of orthologous genes (e.g. human and worm), named "target_original_orthologs.txt", sorted by target gene ID
Example of the tab-delimited data format (target gene ID, original gene name, original gene ID): 

```
ENSG00000000419	dpm-1	WBGene00022044
ENSG00000000938	src-1	WBGene00005077
ENSG00000001167	nfya-1	WBGene00011614
```

Annotation data like these are available from various public annotations sources, e.g. Ensembl, NCBI, or UCSC; Depending on the data source, some reformatting of tables may be possibly needed. For instance, in the two examples above thet data were obtained from [Ensembl Biomart](http://www.ensembl.org/biomart/martview). 

Depending on the data source, some reformatting of tables may be needed; for instance, one may need to make sure human chromosome names are the same in the "target_chip_peaks.bed" and "target_gene_annotation.bed" files above, as conventions may vary between different sources. For Ensembl data, one can add the string prefix "chr" to chromosome names make them compatible with other data; here is an "awk" command to do that:

```
$ awk '{print "chr"$0}' gene_annotation.orig.bed > gene_annotation.updated_chrom_names.bed
```

Other things being equal, it is best to use recent reference genome annotation versions for these data.




#### 7  Ortholog analysis from ChIP-seq data: process data

##### LiftOver

Let us say the target ChIP-seq peak list is in a BED format file named "target_chip_peaks.bed". If the reference genome version of the target ChIP data is different from that of the annotation data, we will need to translate the chromosome coordinates of the ChIP peaks to the that of the annotation reference; for instance we may have to convert the peak BED file from human genome version hg18 to hg38 (the latest human one).

This can be done on the LiftOver program on the UCSC genome browser website:

[http://genome.ucsc.edu/cgi-bin/hgLiftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver)

A user enters the two genomes for the translation (to and from), uploads the BED file and get the a BED file output with the translated coordinates, as well as a file with positions that may have failed translation for some reason (described in that file).
The processing is very fast.

Alternatively, if you have loaded the "ucsc" module on the Odyssey cluster (with "module load ucsc") or have downloaded the [UCSC utilities executables](http://hgdownload.soe.ucsc.edu/admin/exe/), this conversion can be done on the command line.

A LiftOver file, that has the mapping information between the to and from genome versions, has to be downloaded from UCSC which has a separate folder with such files for a given version of a genome. For example, the LiftOver files (with suffixes "over.chain.gz") for human genome version hg18 are under:

[http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver)

To convert from human genome version hg18 to hg38 the LiftOver file needed is: hg18ToHg38.over.chain.gz 

This can be downloaded with:

```
$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
```
The convesrion is done with the UCSC utility program called also "liftOver". For the hg18 to hg38 conversion example the command would be:

```
$ liftOver target_chip_peaks.bed  hg18ToHg38.over.chain.gz  target_chip_peaks.converted.bed  target_chip_peaks.unlifted.bed
```

This is very fast, 1-2 seconds, and returns same results as doing it on the UCSC website. A number of entries (ususally very small) may fail conversion.

The target_chip_peaks.converted.bed file contains the converted peak locations, and entries that failed conversin are in target_chip_peaks.unlifted.bed.


##### Target genes overlaping ChIP peaks

We then sort the converted peak BED file so that chromosome names are sorted alphabetically. This can be done with:

```
$ sort -k1,1d -k2,3n target_chip_peaks.converted.bed | uniq >  target_chip_peaks.converted.sorted.bed
```

The gene annotation file should have been sorted the same way.

Now we have the files we to find target peaks within or up to 1000bp away from target genes, use a bedtools command (bedtools can be loaded with "module load bedtools2" on the cluster). In this command we also add name for each peak in form <chrom_start_end> (eg. chr1_904314_905315) and sort by gene name (first field):

```
$ bedtools window -a target_gene_annotation.bed -b target_chip_peaks.converted.sorted.bed -w 1000 | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$1"_"$6"_"$7"\t"$6"\t"$7}' | sort > target_genes_overlapping_target_chip_peaks.txt
```

Example of format of the resulting file "target_genes_overlapping_target_chip_peaks.txt":

```
ENSG00000000003	chrX	100627109	100639991	chrX_100636214_100637215	100636214	100637215
ENSG00000000005	chrX	100584802	100599885	chrX_100587795_100588796	100587795	100588796
ENSG00000000419	chr20	50934867	50958555	chr20_50957877_50958878	50957877	50958878
```

#####  Orthologs of target genes overlapping ChIP peaks

Can combine our ortholog pairs list with the target genes overlapping peaks to get the original organism (e.g. work ) orthologs corresponding to those target genes with:

```
$ join -t $'\t'  target_original_orthologs.txt   target_genes_overlapping_target_chip_peaks.txt  >  orthologs_of_target_genes_overlapping_ChIP_peaks.txt
```

The output file would like:

```
ENSG00000000419	dpm-1	WBGene00022044	chr20	50934867	50958555	chr20_50957877_50958878	50957877	50958878
ENSG00000000938	src-1	WBGene00005077	chr1	27612064	27635277	chr1_27614277_27615278	27614277	27615278
ENSG00000000938	src-1	WBGene00005077	chr1	27612064	27635277	chr1_27632339_27633340	27632339	27633340
```

Finally, using our "genes_of_interest.txt" list, we can narrow this ortholog list to the orthologs of our genes of interest with:


```
$ for gene in `awk '{print $1}' genes_of_interest.txt`; \
 do \
   awk -v wormgene=$gene '$2==wormgene' orthologs_of_target_genes_overlapping_ChIP_peaks.txt; \
 done | uniq  > gene_of_interest_orthologs_of_target_genes_overlapping_ChIP_peaks.xls
```


#### 8 References

Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2 (2012). Nature Methods 9:357-359.

Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS (2008). Model-based analysis of ChIP-Seq (MACS). Genome Biology 9(9):R137.

Landt SG, Marinov GK, Kundaje A, Kheradpour P, Pauli F, et al. (2012). ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Research 22:1813–1831.

Li Q, Brown J, Huang H, Bickel P (2011). Measuring reproducibility of high-throughput experiments. The Annal of Applied Statistics 5:1752–1779.



