Title: Best Practices for De Novo Transcritome Assembly with Trinity
Date: 2016-08-26 00:00
Author: Adam Freedman
Category: Software
Tags: Next-Gen Sequencing, Transcriptome, Transcriptome Assembly, Trinity
Summary: A best-practices pipeline for de novo transcriptome assembly with Illumina paired-end reads using Trinity


[TRINITY](https://github.com/trinityrnaseq/trinityrnaseq/wiki) is a software package for conducting de novo (as well as the genome-guided version of) transcriptome assembly from RNA-seq data. The Trinity package also includes a number of perl scripts for generating statistics to assess assembly quality, and for wrapping external tools for conducting downstream analyses. 

#### 1 Consult with Informatics Group staff about study design

Not infrequently, Harvard researchers have come to us with assembly problems....AFTER they have devised a study design, picked samples, and conducted sequencing. While in our experience, de novo transcriptome assemblies are far more fragmented than the early performance assessments would suggest, there are aspects of sequencing experiments that can negatively impact assemblies. Please contact us so we can provide advice on designing your transcriptomics study, and help you avoid common pitfalls.  

#### 2 Examine quality metrics for sequencing reads

Examining the base quality distribution, kmer frequencies and adapter contamination by position in the read is an important first step to understanding the underlying quality of your data. For example, an  increase in adapter frequency as one moves along a read is indicative of incomplete removal of adapter sequence during demultiplexing, a rather common occurrence. In addition, the presence of over-represented sequences can be indicative of adapter contamination, rRNA reads, or perhaps other exongenous contamination.

Use [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to examine quality metrics from your reads. An example slurm submission script is as follows: 
                        
    :::bash    
    #!/bin/bash 
    #SBATCH -p serial_requeue       # Partition to submit to 
    #SBATCH -n 1                   # Number of cores 
    #SBATCH -t 0-3:00               # Runtime in days-hours:minutes 
    #SBATCH --mem 2000              # Memory in MB 
    #SBATCH -J FastQC               # job name 
    #SBATCH -o FastQC.%A.out        # File to which standard out will be written 
    #SBATCH -e FastQC.%A.err        # File to which standard err will be written 
    #SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL 
    #SBATCH --mail-user=<PUT YOUR EMAIL ADDRESS HERE>  # Email to which notifications will be sent 

    source new-modules.sh
    module purge
    module load fastqc/0.11.5-fasrc01

    j=`basename $1`

    mkdir -p fastqc_$j

    fastqc --outdir fastqc_$j $1  2>&1 > $j.fastqc.sbatch.out

The above script uses a command line argument, specified by $1, which would be the name of the fastq file. Thus, job submission would look something like:
         
    :::bash
    $ sbatch myfastqcscript.sh mypurpleunicorn_1.fastq  

Outputs will include an html format report, and text files summarizing various quality metric information. Note, in the above script and those below, we use the practice of first loading new-modules to create access to the current module set, then purging any modules that might have gotten loaded from your .bash_profile or .bashrc file to avoid any conflicts, then loading any required modules for the current analysis. 

#### 3  Removing erroneous k-mers from Illumina paired-end reads

Because the current state of the art transcriptome assemblers use a DeBruijn Graph approach that constructs graphs from kmers, erroneous k-mers can adversely impact assemblies. Because rare k-mers are likely due to sequencing errors, correcting reads such that rare k-mers are corrected to a more frequently occurring can improve assemblies. In theory, lowly expressed transcripts may lead to the occurence of biologically real, rare kmers that will get flagged as errors. However, assembly tools will generally do not do a good job of assembling lowly expressed transcripts anyway. Thus, in our opinion the benefits of rooting out errors that will impact the assembly of many transcripts outweigh any adverse effects on reconstruction of lowly expressed transcripts whose assembly will already be compromised by low read coverage. 

We use [rCorrector](https://github.com/mourisl/Rcorrector), a tool specifically designed for RNA-seq data. Besides being a top performer in side-by-side comparisons of kmer-based read error correction tools, includes tags in the fastq output that indicate whether the read has been corrected, or has been detected as containing an error, but is uncorrectable.

First, install a local version of rCorrector. cd into the directory within your home where you install software, and make a clone of rCorrector using git:

   	:::bash
    $ git clone git@github.com:mourisl/Rcorrector.git
    $ make

Then, make a directory in which to run rCorrector, and create symlinks from your fastq files to this directory (to keep your raw data safe). Then, execute an sbatch script:
                
    :::bash
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
     
    perl /path/to/rCorrector/run_rcorrector.pl -t 12 -1 comma-separated list of left (R1) reads -2 comma-separated list of right (R2) reads

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

#### 4  Discard read pairs for which one of the reads is deemed unfixable

Such unfixable reads are often riddled with Ns, or represent other low complexity sequences. Why risk having junk incorporated into your assembly, and why spend compute time on reads that can only hurt the assembly? No good reason, so remove them. We provide a python script for achieving this task. It also strips the "cor" tag from the headers of corrected sequences, as these can cause problems for downstream tools, particularly if you are using data from SRA. Ask not why...it will only give you a headache. 
       
    :::bash
    #!/bin/bash
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

    python FilterUncorrectabledPEfastq.py -1 $1 -2 $2 -o fixed 2>&1 > rmunfixable_$r1.out
     
   # where $1 and $2 are cmd line arguments for names of R1 and R2 fastq files, respectively, and 'fixed' is a prefix added to the filtered output files (which, of course, you can change).

The latest version of [FilterUncorrectablePEfasta.py](https://github.com/harvardinformatics/TranscriptomeAssemblyTools/blob/master/FilterUncorrectabledPEfastq.py) can be found in the TranscriptomeAsemblyTools repository on the Harvard Informatics github repository

#### 5   Trim adapter and low quality bases from fastq files

The Illumina demultiplexing pipeline may incompletely remove adapter sequences, and when the insert sizes for a give read pair lead to overlaps between the sequenced bases, sequencing for one read can extend into the adapter of the other. These bases, as well as low quality bases should be trimmed prior to running Trinity. However, there is evidence that excessive quality trimming (by setting a high base quality threshold) can negatively impact assemblies, [McManes et al. 2014](http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full). Based upon these findings, we currently recommend filtering out bases with qualities below phred 5.  

While Trimmomatic is commonly used for adapter and quality trimming, the adapter composition-by-cycle plots recently added to fastqc have revealed that it does not completely remove adapter sequence, with high rates of adapter contamination remaining in the last bases of a read. In contrast, TrimGalore!, a convenient wrapper for cutadapt, appears to completely remove adapter contamination in reads. Thus, it is our preferred tool. TrimGalore! is not a module on odyssey, but can be downloaded and run out of the box. One merely needs to load the cutadapt module. First, install TrimGalore! into your home directory where you keep software.

    :::bash
    $ wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip
    $ unzip trim_galore_v0.4.1.zip

Then, make and output directory for your trimming results, and submit the trimming job with an sbatch script:
        
    :::bash
    #!/bin/bash
    #SBATCH -J trimgalore
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

    /your/path/to/trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 $1 $2
 
Some of these options can be modified for your data set, e.g. if you are analyzing single end data, you obviously need none of the arguments specifying the handling of paired data, nor do you need to specify R1 and R2 reads! In addition, you may have reasons to change the minimum read length threshold, or the --strigency and -e parameters that pertain to adapter detection and filtering. In our experience, these settings have worked well. If you wish to re-purpose this script for trimming reads that you will then align to a genome, you can opt to be more stringent with respect to the minimum allowed base quality.

**IMPORTANT NOTE:** if any of your downstream applications use bowtie (not bowtie2), it will be necessary to also supply the -t flag to TrimGalore. When a start/end coordinate of one read of an aligned pair is contained within its mate, bowtie will report that alignment as invalid. If you previously used TrimGalore, and then used the older version of the Trinity perl scripts for assessing read support that used bowtie, and did not supply the -t flag, the reported estimate of the percentage of properly mapped pairs will be incorrect, and might substantially under-estimate the support for your assembly.  

Finally, if you have a lot of fastq files that you wish to run separately, in either the fastqc or trimming steps, one should consider writing a loop script that will iterate over files and submit sbatch submissions, rather than manually supplying command line arguments for one pair of reads at a time.

NOTE: For libraries built with Wafergen PrepX directional mRNA library kits on the Apollo robot, we have seen cases where TrimGalore! does not remove adapters in their entirety using TrimGalore's default settings. If the fastqc report for the trimmed reads still indicates an increase in adapter occurrence as one moves along the read, specify the reverse complements of the specific Wafergen adapters and index sequence, such that your TrimGalore! command line for the above submission script looks like this:

        your/path/to/trim_galore --paired --retain_unpaired --phred33 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG -a2 GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT --output_dir trimmed_reads --length 36 -q 5 --stringency 5 -e 0.1 $1 $2

The sequence for -a = A + reverse complement of the index primer (which includes the unique index for the sample). Another way to interpret this is that it is the sequence of the 3' adapter plus the additional bases added to reverse complement the overhanging portion of the index primer. In the above case, the index is ACTTGA, such that expected sequence in the read = AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC + ACTTGA + AATCTCGTATGCCGTCTTCTGCTTG.

The sequenece for -a2 = the reverse complement of the 5' adapter plus the bases added to reverse complement the SR primer. 

See the PrepX workflow documents from Wafergen Biosystems to get more details on the adapter sequences. 

#### 6 Map trimmed reads to a blacklist to remove unwanted (rRNA reads) OPTIONAL

In most cases, researchers will choose poly-A selection over ribo-depletion, either because of interest in coding sequence, or the finding that ribo-depletion can lead to bias in downstream analyses [Lahens et al. 2014, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r86). However, library prepration strategies using poly-A selection will typically not remove all of the rRNA, and we have seen frequencies of reads originating from rRNA post-selection to occasionally exceed 30%. Removing reads originating from rRNA will reduce Odyssey cluster usage, and assembly time. Equally important, you will have a more precise estimate of how many of your reads will actually go towards the assembly of mRNA transcripts.

Our recommendation is to first map your reads to an rRNA database, such as can be downloaded in fasta for from [SILVA](https://www.arb-silva.de/). Then, one can run bowtie2 such as to maximize sensitivity of mapping, meaning you will maximize the number of reads you will consider as originating from rRNA, and thus worthy of being filtered out of your final read set for assembly. Once you build a bowtie2 index for the rRNA fasta database, see instructions at [bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), map your reads:
        
    :::bash
    #!/bin/bash
    #SBATCH -N 1
    #SBATCH -n 12 #Number of cores
    #SBATCH -t 12:00:00  #Runtime in minutes
    #SBATCH -p serial_requeue  #Partition to submit to
    #SBATCH --mem=12000  #Memory per node in MB
    #SBATCH -e silvabt2_%A.e
    #SBATCH -o silvabt2_%A.o
     
    # $1 = full path to your silva database; do not include the fasta suffix (fa,fasta,etc.)
    # $2 = R1 fastq file
    # $3 = R2 fastq file
    # $4 = text file to which to write summary metrics
    # $5 = name of file to which paired reads that map to rRNA database will be written (bowtie2 will and the 1 and 2 suffixes to indicate whether left or right reads
    # $6 = same as $5, but for read pairs for which neither read mapped to the rRNA database
    # $7 = same as $5, but for single ends that mapped to the rRNA database whose mates didn't map
    # $8 = same as $5, but for single ends that did not map to the rRNA database
    
    source new-modules.sh
    module purge
    module load bowtie2/2.2.4-fasrc01
    bowtie2 --very-sensitive-local --phred33  -x $1 -1 $2  -2 $3 --threads 12 --met-file $4 --al-conc-gz $5 --un-conc-gz $6 --al-gz $7 --un-gz $8

The reads you want to keep are those corresponding to the read pairs that did not align to the rRNA database, i.e. specified by $6 after the --un-conc-gz flag.

#### 7 Run fastqc on your processed reads that pass qc and filtering from the above steps

Use the same sbatch submission format from step 2 (above). Ideally, there will be no trend in adapter contamination by cycle, and there will be increased evenness in kmer distributions, GC content, no over-represented sequences, etc. 

#### 8 Remove remaining over-represented sequences OPTIONNAL

Occasionally, over-represented sequeuences will be detected by fastqc even after running through the above steps. One should BLAST them to see what they are, and consider using a script to remove read pairs containing the over-represented sequences. Typically, these will be rRNAs that were not sufficiently represented in the SILVA database (e.g. 5S), or other more common rRNAs that your reads won't map to because it is too divergent from the organisms used to construct the database. We provide a python script,[RemoveFastqcOverrepSequenceReads.py](https://github.com/harvardinformatics/TranscriptomeAssemblyTools/blob/master/RemoveFastqcOverrepSequenceReads.py) to use the sequences flagged by fastqc to remove read pairs where either of the reads contains one of their respective over-represented sequences.

#### 9 Run Trinity
 
We continue to evaluate other de novo transcriptome assemblers, but at present we recommend Trinity as it performs relatively well, uses compute resources efficiently, and has ongoing support from its developers, and the distribution includes scripts for conducting a number of downstream analyses for assembly quality evaluation and expression estimation. 

Settings used for Trinity will depend upon a number of factors, including the sequencing strategy, whether libraries are stranded, and the size of the data set. For data sets with > 200 million reads after the filtering steps above, for computational considerations we recommend conducting in silico normalization. One can perform the normalization as part of the Trinity run, or do it separately using insilico_read_normalization.pl, the perl script provided as part of the Trinity package, found in the util directory. We suggest it may be more efficient to run the normalization separately, especially if a large amount of memory needs to be allocated to the assembly itself. Below is an example script for a Trinity job for a stranded paired-end library, generated using dUTP chemistry (hence the RF flag), without normalization as part of the job submission.
        
    :::bash
    #!/bin/bash 
    #SBATCH -N 1
    #SBATCH -n 32
    #SBATCH -p general                   # may consider running on a bigmem node for large dataset
    #SBATCH -e trinity_%A.err            # File to which STDERR will be written
    #SBATCH -o trinity_%A.out           # File to which STDOUT will be written
    #SBATCH -J trinity_%A               # Job name
    #SBATCH --mem=250000                 # Memory requested
    #SBATCH --time=2-23:00:00              # Runtime in D-HH:MM:SS
    #SBATCH --mail-type=ALL              # Type of email notification- BEGIN,END,FAIL,ALL
    #SBATCH --mail-user=name@harvard.edu # Email to send notifications to

    source new-modules.sh
    module purge
    module load trinity/2.2.0-fasrc01
    # $1 = comma-separated list of R1 files
    # $2 = comma-separated list of R2 files
    # $3 = name of output directory Trinity will create to store results. This must include Trinity in the name, otherwise the job will terminate

    Trinity --seqType fq --SS_lib_type RF --max_memory 225G --min_kmer_cov 1 --CPU 32 --left $1 --right $2 --output $3 

Depending upon the size of the input read data set, you may need some combination of more (or less) time and memory. Another option is to use Trinity's capability to distribute parallelizable steps across a job grid. To do this, you must add the --grid_conf argument and specifying a grid config file, configured as follows:

    # grid type:
    grid=SLURM
    # template for a grid submission
    cmd=sbatch -p general --mem=5500 --time=02:00:00 --account=<YOUR ODYSSEY USER ACCT>
    # number of grid submissions to be maintained at steady state by the Trinity submission system
    max_nodes=500
    # number of commands that are batched into a single grid submission job.
    cmds_per_node=200

You should also use the --grid_node_max_memory 5G argument, which uses a value slightly less than what's in the grid config file. 

Finally, --left and --right are for comma separated lists of R1 and R2 fastq files. An alternative way for specifying a large number of fastq files is to instead use --left_list and --right_list and have the arguments point to txt files that provide the full path names of the R1 and R2 files, respectively, with 1 row per file

