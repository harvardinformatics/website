Title: Metatranscriptome analysis example on Odyssey
Date: 2016-02-08 00:00
Category: Tutorials
Tags: Next-Gen Sequencing, Microbiome, RNA-Seq, Metatranscriptome 
Summary: An example of preprocessing and analyzing microbiome RNA-seq data on the Odyssey cluster

Next generation sequencing (NGS) RNA sequencing (RNA-seq) can be applied in complex microbial ecosystems for metatranscriptome analysis. The application of mRNA enrichment procedures in combination with NGS techniques has established that in depth analysis of transcriptomic landscapes are now feasible in marine and soil microbial communities, [human microbiomes](http://commonfund.nih.gov/hmp/index) from the intestine and many other body parts, and other types of microbial communities. Paired-end sequencing, like that available on the Illumina platform, can be applied to increase the information content and to enhance the identification of expressed gene functions in complex communities.


#### 1  Load the software and start SLURM session

Many software programs are already installed on the cluster as modules and can be loaded using the Lmod system. To activate Lmod use:

	$ source new-modules.sh

Programs can then be loaded with the "module load" command. We start by loading Python; and Mercurial, which is a free, distributed source control management tool; we will need it to install some scripts used in this example workflow:

	$ module load python
	$ module load mercurial 

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), a program for checking the quality of FASTQ read data, is installed as a module on the cluster. It can be loaded as follows:

	$ module load fastqc
	
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) trimming tool for cleaning and filtering Illumina NGS single-end and paired-end reads, can be easily installed locally, by downloading and unzipping the binary executable:

	$ wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip
	$ unzip Trimmomatic-0.35.zip

You can mark the location of the Trimmomatic directory with an environment variable, for easy reference:  
	
	$ export TRIMMOMATIC_HOME=/path/to/trimmomatic/Trimmomatic-0.35

[MetaPhlAn](https://bitbucket.org/biobakery/metaphlan2) is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea, Eukaryotes and Viruses) from metagenomic shotgun sequencing data with species level resolution. It can be installed with (Mercurial) command "hg":

	$ hg clone https://bitbucket.org/biobakery/metaphlan2

Mark the location of the MetaPhlAn directory with an environment variable for easy reference, and export its location in the $PATH environment variable:  

	$ export mpa_dir=/path/to/metaphlan2
	$ export PATH=${mpa_dir}:$PATH


[Krona](https://github.com/marbl/Krona/wiki) is an interactive visualization tool that allows hierarchical data to be explored with zoomable pie charts. 

Download KronaTools and unpack the archive:

	$ wget https://github.com/marbl/Krona/releases/download/v2.6.1/KronaTools-2.6.1.tar
	$ tar xvf KronaTools-2.6.1.tar

Then go to the resulting directory "KronaTools-2.6.1" and issue the following command to install the Krona scripts in your home directory bin ($HOME/bin which is on your $PATH):

	$ ./install.pl --prefix $HOME

We recommend running longer jobs as SLURM jobs on the cluster, or in an interactive SLURM session on Odyssey, for smaller datasets.

When logged in on Odyssey, an interactive session can be launched with a command like the following, in which we are requesting one CPU (parameter "-n") in the interactive partition (parameter "-p"), for a bash session of 150 minutes ("-t") with memory 5000 Mb ("--mem="):

	
	$ srun -p interact -n 1 --pty -t 150 --mem=5000 /bin/bash


Create a separate working folder for this example, for data and results files, running jobs, etc., preferably on the /n/regal file system, which is recommended for SLURM jobs, especially if there is a lot of I/O or large files are written out during jobs.
For example, if you create your working directory under regal and name it "metatrans_example", the full path to the directory could be:

	
	$ export METATRANS_DIR=~/regal/metatrans_example



#### 2  Download RNA-seq metatranscriptome data


If you do not have you own RNA-seq fastq read data, you can download the following for this example.

The RNA-seq example data are from an in vivo sample taken from the anterior nares of a human nose, and are available on GEO at accession [GSE56294](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56294), as GEO Short Read Archive (.sra) datafile with ID "SRR1206356". It can be downloaded with:

	$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX502/SRX502160/SRR1206356/SRR1206356.sra


This should take a few minutes to download "SRR1206356.sra" into the working directory.

To extract fastq files from the GEO data, use NCBI program toolkit. It is already installed on Odyssey and can be loaded with:

	$ source new-modules.sh
	$ module load sratoolkit/0.2.3.4-2.bib-fasrc02

The fastq-dump command compressed fastq files. In this example, this is done with a SLURM job:

	#!/bin/bash
	#SBATCH -n 1
	#SBATCH -N 1
	#SBATCH --mem 1000
	#SBATCH -p serial_requeue
	#SBATCH -o sra.out
	#SBATCH -e sra.err
	#SBATCH -J sra_arr
	#SBATCH -t 50

	source new-modules.sh
	module load sratoolkit/0.2.3.4-2.bib-fasrc02

	METATRANS_DIR=~/regal/metatrans_example

	srafile=SRR1206356

	fastq-dump --defline-seq '@$sn[_$rn]/$ri' --gzip  -O $METATRANS_DIR  $METATRANS_DIR/$srafile.sra

The fastq-dump processing of this should take a few minutes. For small .sra files, an interactive SLURM session would also have worked. For larger input files, use a SLURM job.



#### 3  Preprocess the data

Check quality statistics of reads using FastQC. This is done here using a SLURM job:

	#!/bin/bash   
	#
	#SBATCH -n 1  
	#SBATCH -N 1  
	#SBATCH -t 50  
	#SBATCH -p serial_requeue 
	#SBATCH --mem=3000  
	#SBATCH -o fastqc.out  
	#SBATCH -e fastqc.err 
	#SBATCH -J fastqc

	source new-modules.sh
	module load fastqc

	METATRANS_DIR=~/regal/metatrans_example
	srafile=SRR1206356

	mkdir $METATRANS_DIR/fastqc_out
	fastqc --casava --noextract --nogroup --outdir $METATRANS_DIR/fastqc_out ${srafile}_1.fastq.gz

The FastQC results can be viewed in .html results file which is found in the output folder we created: $METATRANS_DIR/fastqc_out


If there are many reads with low quality bases, trim with Trimmomatic using the following SLURM job:

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
	METATRANS_DIR=~/regal/metatrans_example
	srafile=SRR1206356

	java -jar $TRIMMOMATIC_HOME/trimmomatic-0.35.jar SE -threads 8 -phred33 \
	 $METATRANS_DIR/${srafile}_1.fastq.gz $METATRANS_DIR/${srafile}_1.trimmed.fastq.gz \
	 ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36




#### 4  Analyze the metatranscriptome data

Run MetaPhlAn to profile the microbial composition of the metatranscriptome dataset (in this example as a SLURM job): 

	#!/bin/bash
	#
	#SBATCH -n 8
	#SBATCH -N 1
	#SBATCH --mem=3000
	#SBATCH -p serial_requeue
	#SBATCH -o mpa.out
	#SBATCH -e mpa.err
	#SBATCH -J mpa_job
	#SBATCH -t 100

	source new-modules.sh
	module load bowtie2

	export PATH=/path/to/metaphlan2:$PATH
	mpa_dir=/path/to/metaphlan2

	srafile=SRR1206356

	metaphlan2.py --input_type fastq -o abundance_profile.txt --mpa_pkl ${mpa_dir}/db_v20/mpa_v20_m200.pkl \
	   --bowtie2db ${mpa_dir}/db_v20/mpa_v20_m200 --bowtie2out metagenome.bowtie2.bz2 \
	   --nproc 8 <(gunzip -c ${srafile}_1.trimmed.fastq.gz)

The script produces output file "abundance_profile.txt" (as specified in the metaphlan2.py command parameter "-o").
The MetaPhlAn2 output (microbial abundance tables in tab-delimited format) contains the microbial species (Column 1) and their associated relative abundances (Column 2).

Here are the first few lines of the output, with the taxonomic labels and the percentages for each:

	#SampleID       Metaphlan2_Analysis
	k__Viruses      82.13057
	k__Bacteria     17.86943
	k__Viruses|p__Viruses_noname    82.13057
	k__Bacteria|p__Actinobacteria   11.51396
	k__Bacteria|p__Firmicutes       6.35547
	k__Viruses|p__Viruses_noname|c__Viruses_noname  82.13057
	k__Bacteria|p__Actinobacteria|c__Actinobacteria 11.51396
	k__Bacteria|p__Firmicutes|c__Bacilli    6.10394
	k__Bacteria|p__Firmicutes|c__Clostridia 0.25153
	k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname        82.13057
	k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales      11.51396
	k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales 4.04565
	k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales      2.05829
	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales        0.25153
	k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Potyviridae 34.06371
	k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Parvoviridae        26.06781
	k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Partitiviridae      16.84973


To view the microbial abundance profile graphically, use Krona as follows. First, using a MetPhlAn2 utility, generate an input file abundance_profile.krona.txt for KronaTools:

	$ python $mpa_dir/utils/metaphlan2krona.py -p abundance_profile.txt -k abundance_profile.krona.txt

Then, run the KronaTools command "ktImportText" on this text file to generate an interactive html file ("abundance_profile.krona.html"):

	$ ktImportText -o abundance_profile.krona.html abundance_profile.krona.txt 

A screen capture of the generated html file is shown below.

![abundance_profile.krona]({static}/images/abundance_profile.krona.png)  



#### 5  Metatranscriptome assembly with Trinity

To assemble transcripts from the metatranscriptome reads we can use Trinity. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) is a method for efficient and robust de novo reconstruction of transcriptomes from RNA-seq data. In this example, we are running Trinity for single-end reads ("--single" parameter), 16 CPU cores ("--CPU" parameter), maximum memory of 50GB ("--max_memory" parameter) and minimum contig length 250bp ("--min_contig_length" parameter).
Trinity writes out large intermediate files, so make sure you run the jobs on appropriate filesystems, like /n/regal as was mentioned above.

	#!/bin/bash
	#
	#SBATCH -n 16
	#SBATCH -N 1
	#SBATCH --mem=50000
	#SBATCH -p serial_requeue
	#SBATCH -o trinity.out
	#SBATCH -e trinity.err
	#SBATCH -J trinity_job
	#SBATCH -t 1400
	
	source new-modules.sh
	module load trinity
	
	srafile=SRR1206356

	Trinity --seqType fq --max_memory 50G --single ${srafile}_1.trimmed.fastq.gz --CPU 16 --min_contig_length 250 --output ${srafile}.trinity_results

The results are in output folder with ".trinity_results" suffix; the main output is the Fasta with the assembled transcripts, "Trinity.fasta".

In order to annotate these transcripts, you can use the [Trinotate](https://trinotate.github.io/) tool (from the same group who produced Trinity)<!--: an example workflow for Trinotate can be found in this [document](trinotate-workflow-example-on-odyssey.html)-->. 


#### 6 References

Truong DT, Franzosa EA, Tickle TL, Scholz M, Weingart G, Pasolli E, Tett A, Huttenhower C and Segata N (2015). MetaPhlAn2 for enhanced metagenomic taxonomic profiling. Nature Methods 12, 902–903.

Ondov BD, Bergman NH, and Phillippy AM (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics 12(1):385.

Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A (2011). Full-length transcriptome assembly from RNA-seq data without a reference genome. Nature Biotechnology 29(7):644-52. 



