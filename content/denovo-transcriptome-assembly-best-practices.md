Title: Best Practices for De Novo Transcriptome Assembly with Trinity
Date: 2020-09-29 00:00
Modified: 2021-03-05
Author: Adam Freedman, Nathan Weeks
Category: Tutorials
Tags: Next-Gen Sequencing, Transcriptome, Transcriptome Assembly, Trinity
Summary: A best-practices pipeline for de novo transcriptome assembly with Illumina paired-end reads using Trinity

Despite a steady increase in the availablity of tools and documented pipelines for building transcriptome assemblies, de novo transcriptome assembly from relative short Illumina paired-end reads remains an extremely challenging endeavor. While early genome assemblers used pairwise overlaps between long reads to extend contigs, this approach is unfeasible when dealing with hundreds of millions of reads. Thus, de novo transcriptome asssemblers use DeBruijn graphs, which are constructed and extended based upon kmers, i.e. subsequence of length k found in reads. While this makes the assembly process computationally tractable, it can lead to fragmented assemblies of a large number of contigs that are subsequences of the underlying true transcripts. Some of the factors that lead to this fragmentation are sequencing errors, polymorphism, sequence repeats, and for more lowly expressed transcripts, stochasticicty of read depth that leads to gaps in coverage. Transcriptome assemblers, unlike genome assemblers, must handle the wide range of depth of coverage due to gene expression variation. Our goal in developing a best practice pipeline is to produce most contiguous, error-free and complete transcriptome assemblies given these challenges. While a specific pipeline may produce near optimal results in most scenarios, there may be particular biological and technical factors worth considering that may lead to modest changes to what we present below. As always, it is best to consult with Informatics Group bioinformaticians to discuss these issues.   


[TRINITY](https://github.com/trinityrnaseq/trinityrnaseq/wiki) is a software package for conducting de novo (as well as the genome-guided version of) transcriptome assembly from RNA-seq data. The Trinity package also includes a number of perl scripts for generating statistics to assess assembly quality, and for wrapping external tools for conducting downstream analyses. 

#### 1 Consult with Informatics Group staff about study design

Not infrequently, Harvard researchers have come to us with assembly problems....AFTER they have devised a study design, picked samples, and conducted sequencing. While in our experience, de novo transcriptome assemblies are far more fragmented than the early performance assessments would suggest, there are aspects of sequencing experiments that can negatively impact assemblies. Please contact us so we can provide advice on designing your transcriptomics study, and help you avoid common pitfalls.  

#### 2 Examine quality metrics for sequencing reads

Examining the base quality distribution, kmer frequencies and adapter contamination by position in the read is an important first step to understanding the underlying quality of your data. For example, an  increase in adapter frequency as one moves along a read is indicative of incomplete removal of adapter sequence during demultiplexing, a rather common occurrence. In addition, the presence of over-represented sequences can be indicative of adapter contamination, rRNA reads, or perhaps other exongenous contamination.

Use [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to examine quality metrics from your reads. An example slurm submission script is as follows: 
                        
    :::bash    
    #!/bin/bash 
    #SBATCH -p serial_requeue       # Partition to submit to 
    #SBATCH -N 1
    #SBATCH -n 6                    # Number of cores
    #SBATCH -t 0-3:00               # Runtime in days-hours:minutes 
    #SBATCH --mem 6000              # Memory in MB
    #SBATCH -J FastQC               # job name 
    #SBATCH -o FastQC.%A.out        # File to which standard out will be written 
    #SBATCH -e FastQC.%A.err        # File to which standard err will be written 
    #SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL 
    #SBATCH --mail-user=<PUT YOUR EMAIL ADDRESS HERE>  # Email to which notifications will be sent 

    readonly SINGULARITY_IMAGE='singularity exec --cleanenv /n/singularity_images/informatics/trinityrnaseq/trinityrnaseq.v2.12.0.simg'

    ${SINGULARITY_EXEC} fastqc --threads ${SLURM_CPUS_ON_NODE} --outdir `pwd` $1

The above script uses a command line argument, specified by $1, which would be the name of the fastq file. Thus, job submission would look something like:
         
    :::bash
    $ sbatch myfastqcscript.sh mypurpleunicorn_1.fastq  

Outputs will include an html format report, and text files summarizing various quality metric information. Note, in the above script and those below, we use the practice of first loading new-modules to create access to the current module set, then purging any modules that might have gotten loaded from your .bash_profile or .bashrc file to avoid any conflicts, then loading any required modules for the current analysis. 

#### 3  Removing erroneous k-mers from Illumina paired-end reads

Because the current state of the art transcriptome assemblers use a DeBruijn Graph approach that constructs graphs from kmers, erroneous k-mers can adversely impact assemblies. Because rare k-mers are likely due to sequencing errors, correcting reads such that rare k-mers are corrected to a more frequently occurring can improve assemblies. In theory, lowly expressed transcripts may lead to the occurence of biologically real, rare kmers that will get flagged as errors. However, assembly tools will generally do not do a good job of assembling lowly expressed transcripts anyway. Thus, in our opinion the benefits of rooting out errors that will impact the assembly of many transcripts outweigh any adverse effects on reconstruction of lowly expressed transcripts whose assembly will already be compromised by low read coverage. 

We use [rCorrector](https://github.com/mourisl/Rcorrector), a tool specifically designed for RNA-seq data. Besides being a top performer in side-by-side comparisons of kmer-based read error correction tools, includes tags in the fastq output that indicate whether the read has been corrected, or has been detected as containing an error, but is uncorrectable.

Make a directory in which to run rCorrector, and create symlinks from your fastq files to this directory (to keep your raw data safe). Then, execute an sbatch script:
                
    :::bash
    #!/bin/bash 
    #SBATCH -N 1
    #SBATCH -n 12
    #SBATCH -p general,shared
    #SBATCH -e rcorrect_%A.e
    #SBATCH -o rcorrect_%A.o
    #SBATCH -J rcorrect
    #SBATCH --mem=24000
    #SBATCH --time=36:00:00
    
    module purge
    module load Rcorrector/20180919-fasrc01 
    perl /n/helmod/apps/centos7/Core/Rcorrector/20180919-fasrc01/bin/run_rcorrector.pl -t 12 -1 $1 -2 $2

The $1 and $2 in the script are command line arguments for comma-separated lists of R1 and R2 files for your paired-end reads. Thus, if the script is named rcorrector.sh you would submit the job as follows:

    :::bash
    $ sbatch rcorretor.sh mypurpleunicorn_1.fastq,myyellowunicorn_1.fastq mypurpleunicorn_2.fq,myyellowunicorn_2.fq 

Output fastq files will include "cor" in their names. Corrected reads will have a "cor" suffix in their labels, e.g.

    @ERR1101637.57624042/1 l:165 m:203 h:218 cor
    GTACAACCCTTCCAACCTCCACCGTCTTATATACGAAGCGCCTTGAGTGTGTGTGTGCATGAGCCAAAGGGAATACCG
    +
    <@@D=D2ACDAHB?:<8EGE;B;FE<?;??DBDD))0)8@GDF@C)))5-;5CAE=..7;@DEEC@C;A9?BB=@>@B        

The l,m, and h values in the read headers indicate the kmer counts (across the set of fastq files) for the lowest, median, and highest frequency kmers occuring in the read.

Reads with erroneous kmers that cannot be fixed will also be flagged, e.g.

    @ERR1101637.65/1 l:1 m:45 h:71 unfixable_error
    GGTTCTGTGCCTTTCCCTACCGGAGAGGCGCCCCCAAGACGTATAAGAAGCACCGCTCTTCTCACGCGCATACCAACAACTTCGCAAAGTCTGTGGTCAAC
    +
    @<@DDFBDHHBHHGIIIIIIIIGAG:EGAEGEIIIIBGGHI;C>EAA).;BBDDCC>BD@>>CCABBDD@-2@?::8(2?@@AB<@<@C>CC>A>?944>@

#### 4  Discard read pairs for which one of the reads is deemed unfixable

Such unfixable reads are often riddled with Ns, or represent other low complexity sequences. Why risk having junk incorporated into your assembly, and why spend compute time on reads that can only hurt the assembly? No good reason, so remove them. We provide a python script for achieving this task. It also strips the "cor" tag from the headers of corrected sequences, as these can cause problems for downstream tools, particularly if you are using data from SRA. Ask not why...it will only give you a headache. A python script to remove any read pair where at least one read has an unfixable error, can be obtained from the Harvard Informatics GitHub repository [TranscriptomeAssemblyTools](https://github.com/harvardinformatics/TranscriptomeAssemblyTools). Clone the repository, and run the script in SLURM using the follwoing script: 
       
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

    module purge
    module load python/2.7.8-fasrc01

    python /PATH/TO/FilterUncorrectabledPEfastq.py -1 $1 -2 $2 -s $3
     
The $1 and $2 are cmd line arguments for names of R1 and R2 fastq files, respectively, and, $3 is a sample name to be added to the log that summarized how many reads were removed. 


#### 5   Trim adapter and low quality bases from fastq files

The Illumina demultiplexing pipeline may incompletely remove adapter sequences, and when the insert sizes for a give read pair lead to overlaps between the sequenced bases, sequencing for one read can extend into the adapter of the other. These bases, as well as low quality bases should be trimmed prior to running Trinity. However, there is evidence that excessive quality trimming (by setting a high base quality threshold) can negatively impact assemblies, [McManes et al. 2014](http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full). Based upon these findings, we currently recommend filtering out bases with qualities below phred 5.  

While Trimmomatic is commonly used for adapter and quality trimming, the adapter composition-by-cycle plots recently added to fastqc have revealed that it does not completely remove adapter sequence, with high rates of adapter contamination remaining in the last bases of a read. In contrast, TrimGalore!, a convenient wrapper for cutadapt, appears to completely remove adapter contamination in reads. Thus, it is our preferred tool. TrimGalore! is not a module on odyssey, but can be downloaded and run out of the box. One merely needs to load the cutadapt module. First, install TrimGalore! into your home directory where you keep software.

    :::bash
    $ wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.0.zip
    $ unzip trim_galore_v0.6.0.zip

Then, make and output directory for your trimming results, create symlinks of the Rcorrector-corrected fastq files (with unfixable reads removed), and submit the trimming job with an sbatch script:
        
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
 
    module purge
    module load cutadapt/1.8.1-fasrc01

    # $1 = R1 reads
    # $2 = R2 reads

    /your/path/to/trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 $1 $2
 
Some of these options can be modified for your data set, e.g. if you are analyzing single end data, you obviously need none of the arguments specifying the handling of paired data, nor do you need to specify R1 and R2 reads! In addition, you may have reasons to change the minimum read length threshold, or the --strigency and -e parameters that pertain to adapter detection and filtering. In our experience, these settings have worked well. If you wish to re-purpose this script for trimming reads that you will then align to a genome, you can opt to be more stringent with respect to the minimum allowed base quality.

**Note to bowtie users:** if any of your downstream applications use bowtie (not bowtie2), it will be necessary to also supply the -t flag to TrimGalore. When a start/end coordinate of one read of an aligned pair is contained within its mate, bowtie will report that alignment as invalid. If you previously used TrimGalore, and then used the older version of the Trinity perl scripts for assessing read support that used bowtie, and did not supply the -t flag, the reported estimate of the percentage of properly mapped pairs will be incorrect, and might substantially under-estimate the support for your assembly.  

**Libraries build with Wafergen PrepX mRNA kit on Apollo robot:** For libraries built with Wafergen PrepX directional mRNA library kits on the Apollo robot, we have seen cases where TrimGalore! does not remove adapters in their entirety using TrimGalore's default settings. If the fastqc report for the trimmed reads still indicates an increase in adapter occurrence as one moves along the read, specify the reverse complements of the specific Wafergen adapters and index sequence, such that your TrimGalore! command line for the above submission script looks like this:

        your/path/to/trim_galore --paired --retain_unpaired --phred33 -a AAGATCGGAAGAGC  -a2 GATCGTCGGACTGTAGAA --output_dir trimmed_reads --length 36 -q 5 --stringency 5 -e 0.1 $1 $2

The sequence supplied for -a is A + a subsequences of the reverse complement of the 3' index primer (including, when used, a barcode index). For -a2, it is a subsequences of the reverse complement of the 5' index primer. While, in principle, one could supply the full length adapter sequences (which would include the reverse complements of overhangs), in practice we have found that using the shorter sequences increases sensitivity, and minimizes the chance of retaining adapter sequence in trimmed reads. 

See the PrepX workflow documents from Wafergen Biosystems to get more details on the adapter sequences. 

Finally, if you have a lot of fastq files that you wish to run separately, in either the fastqc or trimming steps, one should consider writing a loop script that will iterate over files and submit sbatch submissions, rather than manually supplying command line arguments for one pair of reads at a time.

#### 6 Map trimmed reads to a blacklist to remove unwanted (rRNA reads) -- OPTIONAL

In most cases, researchers will choose poly-A selection over ribo-depletion, either because of interest in coding sequence, or the finding that ribo-depletion can lead to bias in downstream analyses [Lahens et al. 2014, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r86). However, library prepration strategies using poly-A selection will typically not remove all of the rRNA, and we have seen frequencies of reads originating from rRNA post-selection to occasionally exceed 30%. Removing reads originating from rRNA will reduce Cannon cluster usage, and assembly time. Equally important, you will have a more precise estimate of how many of your reads will actually go towards the assembly of mRNA transcripts.

Our recommendation is to first map your reads to an rRNA database, such as can be downloaded in fasta for from [SILVA](https://www.arb-silva.de/). Then, one can run bowtie2 such as to maximize sensitivity of mapping, meaning you will maximize the number of reads you will consider as originating from rRNA, and thus worthy of being filtered out of your final read set for assembly. From SILVA, we download the SSUParc and LSUParc fasta files, concatenating them, and replacing U characters with T, as our sequence reads are in DNA space. Then, for the concatenated fasta file, you build a bowtie2 index for the rRNA fasta database, see instructions at [bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), and map your reads:
        
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
    # $4 = sample_id (no spaces)
    
    singularity exec --cleanenv /n/singularity_images/informatics/trinityrnaseq/trinityrnaseq.v2.12.0.simg bowtie2 --quiet --very-sensitive-local --phred33  -x $1 -1 $2 -2 $3 --threads 12 --met-file ${4}_bowtie2_metrics.txt --al-conc-gz blacklist_paired_aligned_${4}.fq.gz --un-conc-gz blacklist_paired_unaligned_${4}.fq.gz  --al-gz blacklist_unpaired_aligned_${4}.fq.gz --un-gz blacklist_unpaired_unaligned_${4}.fq.gz

Both an R1 and R2 will get written for each of the switches.

The reads you want to keep are those corresponding to the read pairs for which neither read mapped to the rRNA database,i.e. "paired_unaligned" reads, specified  after the --un-conc-gz flag.

**If your RNA-seq libraries are built with a stranded protocol**, you should use the relevant bowtie2 setting:

* for dUTP-based libraries (Illumina TruSeq,NEBNext Ultra Directional), which are "RF" in Trinity parlance,  use bowtie2 flag --nofw
* for ligation-stranded protocols, i.e. Wafergen PrepX directional mRNA libraries built on the Apollo robot, which are "FR" in Trinity parlance, use bowtie2 flag --norc
* if in doubt about which stranding protocol is used in your library prep kit, consult the user manual or contact the manufacturer's technical support
* of course, if your kit is not a strand-specific protocol, ignore all of this!

NOTE: One can build a blacklist database out of any non-target sequences, e.g. parasites,bacteria, other potential sources of exogenous RNA. One can also augment an rRNA database with known rRNA sequences for your taxon or a closely related one, if they do not appear in SILVA.


#### 7 Run fastqc on your processed reads that pass qc and filtering from the above steps

Use the same sbatch submission format from step 2 (above). Ideally, there will be no trend in adapter contamination by cycle, and there will be increased evenness in kmer distributions, GC content, no over-represented sequences, etc. 

#### 8 Remove remaining over-represented sequences -- OPTIONAL

Occasionally, over-represented sequeuences will be detected by fastqc even after running through the above steps. One should BLAST them to see what they are, and consider using a script to remove read pairs containing the over-represented sequences. Typically, these will be rRNAs that were not sufficiently represented in the SILVA database (e.g. 5S), or other more common rRNAs that your reads won't map to because it is too divergent from the organisms used to construct the database. We provide a python script,[RemoveFastqcOverrepSequenceReads.py](https://github.com/harvardinformatics/TranscriptomeAssemblyTools/blob/master/RemoveFastqcOverrepSequenceReads.py) to use the sequences flagged by fastqc to remove read pairs where either of the reads contains one of their respective over-represented sequences.

#### 9 Run Trinity
 
We continue to evaluate other de novo transcriptome assemblers, but at present we recommend Trinity as it performs relatively well, uses compute resources efficiently, and has ongoing support from its developers, and the distribution includes scripts for conducting a number of downstream analyses for assembly quality evaluation and expression estimation. 

Settings used for Trinity will depend upon a number of factors, including the sequencing strategy, whether libraries are stranded, and the size of the data set.

**Normalization.** In the latest version of Trinity, [in silico normalization](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Insilico-Normalization) is done by default. For data sets with > 200 million reads after the filtering steps above, for computational considerations we recommend using the default mode, and allowing Trinity to normalize reads. If you wish to turn off normalization, then include the **--no_normalize_reads** flag in your Trinity command line.


**Running Trinity inside a Singuarlity container image.**
Rebuilding a new Trinity module with each update is increasingly complicated, as functionality gets added along with additional dependencies. Thus, our preferred mode of running Trinity (as well as the [preferred mode for the Trinity developers](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker#running-trinity-using-singularity)) is to do so inside a [Singularity](https://docs.rc.fas.harvard.edu/kb/singularity-on-the-cluster/) container image, which operates like a light-weight virtual machine, within which software dependencies are conveniently bundled, and relevant environment variables are properly set.
For the convenience of users of the Harvard FAS Cannon Cluster, we provide images from the [Trinity Singularity Image Archive](https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/) at /n/singularity_images/informatics/trinityrnaseq/.

---

*NOTE*: The directions in this guide refer to Trinity v2.12.0.
To see all supported versions of Trinity available, use the FAS RC [module-query](https://docs.rc.fas.harvard.edu/kb/modules-intro/) command: `module-query trinityrnaseq`

---


Running Trinity via Singularity involves two steps. First we run Trinity as a SLURM job. Below is an example script for a Trinity job (with normalization):

> **The current recommended file system on Cannon from which to run Trinity jobs is holyscratch01, as it is colocated in the same data center as the Cannon compute nodes (Holyoke).
> Do not submit Trinity jobs from a boslfs / boslfs02 file system (located in Boston).
> Input FASTQ/FASTA files *may* reside on boslfs/boslfs02 as long as the `trinity_out_dir` (created by default in the directory from which the job script is sumitted; see the TRINITY_OUT_DIR variable below) is located on holyscratch01.**
     
    :::bash
    #!/bin/sh
    # use an entire node in the specified Slurm partition
    #SBATCH --nodes=1
    #SBATCH --mem=0
    #SBATCH --exclusive
    # Adjust wall time limit as appropriate for partition. If the job is cancelled
    # due to time limit exceeded, the job script can be resubmitted, and Trinity
    # will resume execution after the last completed.
    #SBATCH --time=72:00:00
    # Resubmit to the "bigmem" partition if inchworm std::bad_alloc error occurs.
    # The job will then stop after inchworm completes, and can be resubmitted to the
    # "shared" or "test" partition, which are preferred for faster / more reliable
    # execution and sooner job start time.
    #SBATCH --partition=shared

    set -o nounset -o errexit -o xtrace

    ########################################
    # parameters
    ########################################
    readonly SINGULARITY_IMAGE=/n/singularity_images/informatics/trinityrnaseq/trinityrnaseq.v2.12.0.simg
    readonly TRINITY_OUT_DIR=trinity_out_dir
    # To see all options:
    #    singularity exec --cleanenv ${SINGULARITY_IMAGE}  Trinity --show_full_usage_info
    readonly TRINITY_OPTIONS="--output ${TRINITY_OUT_DIR} --max_memory $((8*$(ulimit -m)/(1024**2)/10))G --CPU ${SLURM_CPUS_ON_NODE} $@"

    ########################################
    # ... don't modify below here ...
    if [ ! -s "${TRINITY_OUT_DIR}/read_partitions.img" ]
    then
      mkdir -p "${TRINITY_OUT_DIR}"
      readonly tmpdir=$(mktemp -d)
      mkdir -m 777 -p ${tmpdir}/upper ${tmpdir}/work
      truncate -s 2T "${TRINITY_OUT_DIR}/read_partitions.img"
      singularity exec --cleanenv ${SINGULARITY_IMAGE} mkfs.ext3 -d "${tmpdir}" "${TRINITY_OUT_DIR}/read_partitions.img"
      singularity exec --cleanenv --overlay ${TRINITY_OUT_DIR}/read_partitions.img ${SINGULARITY_IMAGE} mkdir /read_partitions
      ln -sf /read_partitions ${TRINITY_OUT_DIR}/read_partitions
      rm -rf "${tmpdir}"
    fi

    # if on a bigmem node, stop after inchworm
    case ${SLURM_JOB_PARTITION} in
      bigmem) no_run_chrysalis='--no_run_chrysalis' ;;
           *) no_run_chrysalis='' ;;
    esac

    srun -n 1 env time -v singularity exec \
                               --cleanenv \
                               --no-home \
                               --overlay ${TRINITY_OUT_DIR}/read_partitions.img \
                               "${SINGULARITY_IMAGE}" \
      Trinity ${TRINITY_OPTIONS} ${no_run_chrysalis}



This script automatically sets the Trinity `--max_memory` and `--CPU` options based on hardware characteristics of the compute node the job is run on.

If this script is saved to a file called _trinity.sh_, we would submit the Slurm job like this:

    :::bash
    sbatch trinity.sh --seqType fq --left {comma-separated R1 fastq files} --right {comma-separated R2 fastq files}

Alternatively, instead of appending options to the `sbatch` command, the **TRINITY_OPTIONS** string in the trinity.sh script could be edited to reflect particular desired features, e.g.:

* Turning off normalization (`--no_normalize_reads`)
* For directional libraries, `--SS_lib_type` should be set to FR or RF for ligation-stranded and dUTP-based library construction, respectively.
* An alternative way for specifying a large number of fastq files is to instead use `--left_list` and `--right_list` and have the arguments point to txt files that provide the full path names of the R1 and R2 files, respectively, with 1 row per file

Then the job script would be submitted without arguments:

    :::bash
    sbatch trinity.sh


Once the Trinity run has successfully completed, one will need to inspect the results, which are written (by default) to _trinity_out_dir/Trinity.fasta_.

The _trinity_out_dir/read_partitions_ directory can contain hundreds of thousands or millions of files from Chrysalis-generated components (de Brujin graphs & partitioned sequence reads) that are input for Butterfly, which generates additional files containing full-length transcripts.
For I/O efficiency, the contents of the _trinity_out_dir/read_partitions_ directory are written to an [overlay image file](https://sylabs.io/guides/3.5/user-guide/persistent_overlays.html#file-system-image-overlay) (trinity_out_dir/read_partitions.img).
The read_partitions.img image file is created as a [sparse file](https://en.wikipedia.org/wiki/Sparse_file), so it initially consumes far less disk space (as indicated by the `du` command, which reports disk usage) than its maximum size (2T):

    :::bash
    $ ls -lh read_partitions.img
    -rw-rw----+ 1 user group 2.0T Dec  5 16:22 read_partitions.img
    $ du -h read_partitions.img
    34G     read_partitions.img

#### 10-1 Assessing assembly quality step 1: basic alignment summary metrics

Metrics such as N50 should never, by themselves, be treated as good indicators of assembly quality. An obvious, if extreme, example, is that if you incorrectly assembly all of your reads into one gigantic contig, your N50 will be very large. However, extremely short N50s, such that they represent a fraction of the expected size of the fragments in your library might indicate other problems. Similarly, the number of "transcripts" and "genes" in your Trinity assembly do not provide any absolute metric of quality. However, the more "genes" Trinity assembles--particularly if they are only a few hundred bases long--the more likely your contigs represent subsequences of actual genes. These caveats aside, you can easily generate N50 statistics, and counts of the number of Trinity contigs in an [interactive job](https://docs.rc.fas.harvard.edu/kb/running-jobs/#Interactive_jobs_and_srun) using the [TrinityStats.pl](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats#the-gene-contig-nx-statistic) script that comes with Trinity.
    
    :::bash
    srun --pty -p shared -t 00:20:00 --mem 500 /bin/bash
    singularity exec --cleanenv /n/singularity_images/informatics/trinityrnaseq/trinityrnaseq.v2.12.0.simg sh -c '$TRINITY_HOME/util/TrinityStats.pl Trinity.fasta' > Trinity_assembly.metrics

#### 10-2 Assesing assembly quality step 2: quantify read support for the assembly

As explained in the [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly) documentation, assembled transcripts may not represent the full complement of paired-end reads. This will occur because, for very short contigs, only one read from paired-end read will align to it. Simply mapping your reads with bowtie (or your aligner of choice) to the transcripts will not shed any insight into this phenomenon as only properly mapped read pairs will be reported. To evalute read support for the assembly is a three step process. First, you build a bowtie2 index for your assembly.


    :::bash
    #!/bin/bash
    #SBATCH -N 1
    #SBATCH -n 4 #Number of cores
    #SBATCH -t 0:00:00  #Runtime in minutes
    #SBATCH -p serial_requeue,shared  #Partition to submit to
    #SBATCH --mem=8000  #Memory per node in MB
    #SBATCH -e bt2build.e
    #SBATCH -o b2build.o
     
    # $1 = your assembly fasta
    assembly_prefix=$(basename $1 |sed 's/.fasta//g')
    
    readonly SINGULARITY_EXEC='singularity exec --cleanenv /n/singularity_images/informatics/trinityrnaseq/trinityrnaseq.v2.12.0.simg'

    ${SINGULARITY_EXEC} bowtie2-build -–threads 4 $1 $assembly_prefix

Next, you map your reads and calculate alignment statistics.


    :::bash
    #!/bin/bash
    #SBATCH -N 1
    #SBATCH -n 16 #Number of cores
    #SBATCH -t 12:00:00  #Runtime in minutes
    #SBATCH -p serial_requeue,shared  #Partition to submit to
    #SBATCH --mem=16000  #Memory per node in MB

    readonly SINGULARITY_EXEC='singularity exec --cleanenv /n/singularity_images/informatics/trinityrnaseq/trinityrnaseq.v2.12.0.simg'

    # $1 name of your assembly (without the .fasta suffix)
    # $2 comma separated list of left read file names
    # $3 comma separated list of right read file names
    
    ${SINGULARITY_EXEC} bowtie2 -p 10 -q --no-unal -k 20 -x $1 -1 $2 -2 $3  2>align_stats.txt| ${SINGULARITY_EXEC} samtools view -@10 -Sb -o bowtie2.bam

The align_stats.txt file will provide info on the percentage of read pairs that mapped concordantly,as well as an overall alignment rate.

**NOTE: you can use the appropriate flag for stranded libraries (see information in section 6 above re: the appropriate flags to use. Also, depending upon the size of your read data set, you may want to specify more or less time for this job. If you need > 24 hours, be sure to specify the shared partition.


#### 10-3 Assesing assembly quality step 3: quantifying completeness

Another metric of assembly quality is evaluating the extent to which it recovers single copy orthologs that are present across higher taxonomic groupings. While in the absence of knowing which transcripts are truly expressed in a sample it is difficult to determine an absolute expectation for recovery of these orthologs, clearly, high numbers of such genes classified as missing in an assembly should be considered a potential red flag. Furthermore, asssemblies based upon the same read data can be evaluated with respect to the numbers of genes that are complete, fragmented, or missing from the assembly.

To assess completeness, we use [BUSCO](http://busco.ezlab.org/). BUSCO requires that you specify a BUSCO (Benchmarking Universal Single-Copy Orthologs) data set. Data sets are located on the Cannon cluster in /n/holyscratch01/external_repos/INFORMATICS/BUSCO/. Contact us if the database you need is not currently in that directory. BUSCO wraps HMMER, and uses Hidden Markov Model profiles to determine whether assembly contigs are orthologus with a particular BUSCO dataset entry. We recommend you build a conda environment that contains BUSCO, and then running it using that environment. First, you build the environment by launching an [interactive job](https://docs.rc.fas.harvard.edu/kb/running-jobs/#Interactive_jobs_and_srun) (e.g. srun --pty -p shared -t 01:00:00 --mem 1000 /bin/bash), then do:

    :::
    module load python
    conda create -n busco -c bioconda busco

Once the environment has been build, you simply activate it in your BUSCO job script:

    :::
    #!/bin/bash 
    #SBATCH -N 1
    #SBATCH -n 16
    #SBATCH -p serial_requeue,shared
    #SBATCH -e BUSCO.err           # File to which STDERR will be written
    #SBATCH -o BUSCO.out         # File to which STDOUT will be written
    #SBATCH -J BUSCO                # Job name
    #SBATCH --mem=6000                     # Memory requested
    #SBATCH --time=04:00:00                # Runtime in HH:MM:SS
    #SBATCH --mail-type=ALL                # Type of email notification- BEGIN,END,FAIL,ALL
    #SBATCH --mail-user=Your.Email.Address # Email to send notifications to
    
    module load python
    source activate busco

    # $1 input fasta file (your assembly, e.g. Trinty.fasta)
    # $2 output directory (BUSCO will prepend run_)
    # $3 lineage directory see lineages to choose from at: http://busco.ezlab.org/

    run_BUSCO.py -c 16 -o $2 -in $1 -l $3 -m transcriptome

    Example submission: sbatch BUSCO.sh Trinity.fasta trinity_BUSCO /n/holyscratch01/external_repos/INFORMATICS/BUSCO/eukaryota_odb9


Besides directories containing HMMER output, translated proteins, and tables of BUSCO hits (and missing BUSCOs), BUSCO outputs a useful short summary table, and example of which is below:

    BUSCO was run in mode: trans

    Summarized benchmarks in BUSCO notation:
	C:70%[D:24%],F:4.2%,M:25%,n:3023

    Representing:
        1387	Complete Single-copy BUSCOs
        750	Complete Duplicated BUSCOs
        128	Fragmented BUSCOs
        758	Missing BUSCOs
        3023	Total BUSCO groups searched

Because multiple transcripts can hit a BUSCO, this can lead to classification of BUSCOs as being "complete duplicated." More informative is the total number of complete BUSCOs (both single and duplicated). 

 
