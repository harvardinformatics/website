Title: De Novo Transcritome Assembly with Trinity: Odyssey
Date: 2016-08-01 00:00
Category: Software
Tags: Next-Gen Sequencing, Transcriptome, Transcriptome Assembly, Trinity
Summary: A best-practices pipeline for de novo transcriptome assembly with Illumina paired-end reads using Trinity


[TRINITY](https://github.com/trinityrnaseq/trinityrnaseq/wiki) is a software package for conducting de novo (as well as the genome-guided version of) transcriptome assembly from RNA-seq data. The Trinity package also includes a number of perl scripts for generating statistics to assess assembly quality, and for wrapping external tools for conducting downstream analyses. 



#### 1  Removing erroneous k-mers from Illumina paired-end reads

Because Trinity uses a DeBruijn Graph approach, constructing graphs from kmers, erroneous k-mers can adversely impact assemblies. Because rare k-mers are likely due to sequencing errors, correcting reads such that rare k-mers are corrected to a more frequently occurring can improve assemblies. We use rCorrector(https://github.com/mourisl/Rcorrector), which besides being a top performer in side-by-side comparisons, includes tags in the fastq output that indicate whether the read has been corrected, or has been detected as containing an error, but is uncorrectable.

First, install a local version of rCorrector. cd into the directory within your home where you install software, and make a clone of rCorrector using git.
       	:::bash
	$ git clone git@github.com:mourisl/Rcorrector.git
        $ make

Then, make a directory in which to run rCorrector, and create symlinks from your fastq files to this directory (to keep your raw data safe). Then, execute an sbatch script:

        #!/bin/bash 
        #SBATCH -N 1
        #SBATCH -n 12
        #SBATCH -p general
        #SBATCH -e rcorrect_%A.e
        #SBATCH -o rcorrect_%A.o
        #SBATCH -J rcorrect
        #SBATCH --mem=24000
        #SBATCH --time=72:00:00
        
        source new-modules.sh
        module purge
        module load perl/5.10.1-fasrc04
        module load perl-modules/5.10.1-fasrc04
         
        perl /path/to/rCorrector/run_rcorrector.pl -t 12 -1 comma-separated list of left (R1) reads -2 comma-separated list of right (R2)reads

Output fastq files will include "cor" in their names. Corrected reads will have a "cor" suffix in their labels, e.g.
        
        @ILLUMINA-D00365:333:HALEVADXX:1:1101:3092:2052 1:N:0:CAGATC cor
        GTTGAGATGTGTAGGAAGGCAAGTTGGGGTTGATTAAGTCCAATTGTAACTATTATTAGGCCAAGTTGACTTGATGTCGAGAAGGCAATAATTTTTTTAATGTCGTTTTGTGTCAGGGCGCAAATTGCTGTAAATGCAGTGGTTATGGCC
        +
        B?BDDD?DHHDFHIIGIIIIIEDFDHHII?FGIIIIHIFIIGGIIHFHIJIIIJIJJIJJGIJJJHGHGGFHHHHHFFFDDDCCBDDDDDDEEEDDDDDDDDCDB@BBBD?CCDCD9?@BDD@5@CCACCDCD@BCC.9<BB<BACCCBD

Reads with erroneous kmers that cannot be fixed will also be flagged, e.g.

        @ILLUMINA-D00365:333:HALEVADXX:1:1101:1439:2208 1:N:0:CAGATC unfixable_error
        GGAGGCTGCTGCTTGCGTTAGAAAGTACTTGGTGGCGGCTTCTGTGGATCGTGGGTGGTGGGTTGTTGAAATAATTGGAATAATTGATAGGGTANTAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATC
        +
        BC@FFFDFGHGHHJJJJJJJJBHIICGIIJJIGHIJJJJIJJJJIIJHHHHEBFD9AB9?BD8@BDDDDDDDDDEEDDCDCCCDDD@DDDDD9?########################################################

#### 2  Discard read pairs for which one of the reads is deemed unfixable

Such unfixable reads are often riddled with Ns, or represent other low complexity sequences. Why risk having junk incorporated into your assembly, and why spend compute time on reads that can only hurt the assembly? No good reason, so remove them. We provide a python script for achieving this task. It also strips the "cor" tag from the headers of corrected sequences, as these can cause problems for downstream tools, particularly if you are using data from SRA. Ask not why...it will only give you a headache. 
       
        #!/bin/sh
        #SBATCH -J rmunfix_$1                # Job name
        #SBATCH -n 1                         # Use 1 core for the job
        #SBATCH -t 03:00:00                  # Runtime in HH:MM:SS
        #SBATCH -p serial_requeue            # Partition to submit to
        #SBATCH --mem=2000                   # Memory per node in MB
        #SBATCH -o rmunfix_%A.o              # File to which STDOUT will be written
        #SBATCH -e rmunfix_%A.e              # File to which STDERR will be written
        #SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
        #SBATCH --mail-user=name@harvard.edu # Email to which notifications will be sent

        r1=`basename $1`
        r2=`basename $2`

        module purge
        module load python/2.7.8-fasrc01

        python /n/regal/informatics/workshops/trinity/python_code/FilterUncorrectabledPEfastq.py $1 $2 2>&1 > rmunfixable_$r1.out
         
       # where $1 and $2 are cmd line arguments for names of R1 and R2 fastq files, respectively.

#### 3   Trim adapter and low quality bases from fastq files

The Illumina demultiplexing pipeline may incompletely remove adapter sequences, and when the insert sizes for a give read pair lead to overlaps between the sequenced bases, sequencing for one read can extend into the adapter of the other. These bases, as well as low quality bases should be trimmed prior to running Trinity. However, there is evidence that excessive quality trimming (by setting a high base quality threshold) can negatively impact assemblies (CITE McMANES). Based upon these findings, we recommend filtering out bases with qualities below phred 5. 

While Trimmomatic is commonly used for adapter and quality trimming, the adapter composition-by-cycle plots recently added to fastqc have revealed that it does not completely remove adapter sequence, with high rates of adapter contamination remaining in the last bases of a read. In contrast, TrimGalore!, a convenient wrapper for cutadapt, appears to completely remove adapter contamination in reads. Thus, it is our preferred tool. TrimGalore! is not a module on odyssey, but can be downloaded and run out of the box. One merely needs to load the cutadapt module. First, install TrimGalore! into your home directory where you keep software.

       :::bash
       $ wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip
       $ unzip trim_galore_v0.4.1.zip

Then, make and output directory for your trimming results, and submit the trimming job with an sbatch script:

      #!/bin/bash
      #SBATCH -J trimglaore
      #SBATCH -n 1                     # Use 1 cores for the job
      #SBATCH -t 0-6:00                 # Runtime in D-HH:MM
      #SBATCH -p serial_requeue         # Partition to submit to
      #SBATCH --mem=3000               # Memory pool for all cores (see also --mem-per-cpu)
      #SBATCH -o trimgalore_PE.%A.out  # File to which STDOUT will be written
      #SBATCH -e trimgalore_PE.%A.err  # File to which STDERR will be written
      #SBATCH --mail-type=ALL           # Type of email notification- BEGIN,END,FAIL,ALL
      #SBATCH --mail-user=name@harvard.edu # Email to send notifications to
     
      source new-modules.sh
      module purge
      module load cutadapt/1.8.1-fasrc01

      # $1 = R1 reads
      # $2 = R2 reads

      /your/path/to/trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 -q 5 --stringency 5 -e 0.1 $1 $2
 
 Some of these options can be modified for your data set, e.g. if you are analyzing single end data, you obviously need none of the arguments specifying the handling of paired data, nor do you need to specify R1 and R2 reads! In addition, you may have reasons to change the minimum read length threshold, or the --strigency and -e parameters that pertain to adapter detection and filtering. In our experience, these settings have worked well. If you wish to re-purpose this script for trimming reads that you will then align to a genome, you can opt to be more stringent with respect to the minimum allowed base quality.

Finally, if you have a lot of fastq files that you wish to run separately, in either the fastqc or trimming steps, one should consider writing a loop script that will iterate over files and submit sbatch submissions, rather than manually supplying command line arguments for one pair of reads at a time.

#### Running Trinity 

