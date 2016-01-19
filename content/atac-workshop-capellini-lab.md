Title: ATAC Workshop - Capellini Lab
Date: 2015-01-01
Category: Tutorials
Author: Michele Clamp
Tags: Chip-Seq
Summary: This tutorial provides useful tips for Chip-Seq experiment analysis.

## Notes
The sbatch files will appear on the filesystem at

    /n/regal/informatics/workshops/Capellini/Slurm

The output files will appear on the filesystem at

    /n/regal/informatics/workshops/Capellini/Output

## Prerequisites

*   Ability to login to odyssey (ssh, openauth)
*   Knowledge of basic command line commands (ls, cd, cat, less, nano, mkdir, |, >)
*   Knowledge of basic script loops
*   Knowledge of slurm batch scripts
*   Ability to transfer files to/from the cluster

## ATAC-Seq Summary

Outline of ATAC-Seq Sequencing Analysis Process

*   Design experiment and sequence
*   Pre-process the reads (trimmomatic - we're going to skip this step)
*   Align the reads to a reference genome  (bowtie2)
*   Call peaks on the aligned reads using condition and control alignments (macs)
*   Investigate the location of the peaks for enriched annotation and go terms (homer)
*   Extra: Investigate the content of the peaks for binding motifs (homer)
*   Extra: Search the peaks for a specific motif (homer)

## Experimental Design

When starting to think about an ATAC-Seq sequencing experiment several things need to be considered.  These are :

*   number of replicates
*   suitable controls
*   number and type of reads

Guidelines on these questions are addressed in the following publications :

[ENCODE guidelines](http://www.ncbi.nlm.nih.gov/pubmed/22955991) Landt et al 2010 

### Sequencing depth
The number of reads needed depends on the type and frequency of peak expected (sharp/broad,  common/rare,   long motif/short motif) and also the size of the genome. Rough guidelines are as follows:

#### Mammalian cells

Sharp peaks (TFs) -  10 million uniquely mapped (with 80% distinct reads) per replicate Broad peaks (chromatin mods e.g. H3H27me3 M3K36me3) -  12-20 million uniquely mapped (with 80% distinct reads) per replicate  

#### Flies/Worms

Sharp peaks (TFs)</div>

* 2 million uniquely mapped (with 80% distinct reads) per replicate

Broad peaks (chromatin mods e.g. H3H27me3 M3K36me3)

* 5-10 million uniquely mapped (with 80% distinct reads) per replicate


### Controls

These are essential to eliminate false positives as there is always some signal on open chromatin.

### Replicates

The good news here is that, unlike RNA-Seq, more than 2 replicates does not significantly increase the number of targets.

## Workshop dataset

There are 8 files for 4 samples (one file for read1 and one for read2).   We desire to know the difference in  chromatin availability between the distal and proxminal ends of the mouse femur. The files are :   

    Proximal-Femur-1_S1.R1.fastq 
    Proximal-Femur-1_S1.R2.fastq
    Proximal-Femur-3_S3.R1.fastq
    Proximal-Femur-3_S3.R2.fastq
    Distal-Femur-2_S2.R1.fastq
    Distal-Femur-2_S2.R2.fastq
    Distal-Femur-4_S4.R1.fastq
    Distal-Femur-4_S4.R2.fastq 

The raw fastq files are available on the cluster in

    /n/informatics/workshops/Capellini/Data


## Logging on and preparing the data

1.  Log onto one of the cluster login nodes login.rc.fas.harvard.edu
<div class="highlight"><pre>ssh login.rc.fas.harvard.edu</pre></div>

1.  Change directory to our workshop directory `/n/regal/informatics/workshops/Capellini`
<div class="highlight"><pre>cd /n/regal/informatics/workshops/Capellini</pre></div>

1.  Make your own working directory under the Users directory.  You can use your lab storage but regal should be the fastest disk we have.
<div class="highlight"><pre>mkdir Users/mclamp</pre></div>

### Prepare the datasets for alignment. 

#### Symlink the fastq files 
We're going to symlink the fastq files rather than copy them or refer to them in another directory.  That way when you come back to your data in 3 months time you can look in this directory and know which files the alignments came from. 

We are using a shell script loop to save typing out all the filenames.

    cd Users/mclamp
    for i in ../../Data/*.fastq ; do
      ln -s $i .
    done

<div class="exercises">
Exercise<br/>
In your user directory use a for loop to make symlinks to the fastq files.  Use the <code>ls  -l</code> command to check you have links to files and use the head command on all of them to check all the links are correct. (`head *.fastq` will print out the first 10 lines of every file ending fastq) 
</div>

#### The module system.
Odyssey contains a lot of software and by default most bioinformatics tools won't run out of the box if you type them on the command line.    To 'load' them up you need to use the module system.  For instance for bowtie2 and samtools we need to use

<pre>module load centos6/bowtie2-2.2.1
module load centos6/samtools-0.1.19</pre>

If you type in these commands and then try typing bowtie2 -h or samtools you should get some help text. If you want to find  a module use

<pre>module avail 2>&1 |grep -i <searchstring></pre>

For instance to search for samtools

<pre>module avail 2>&1 |grep -i samtools</pre>

Other useful module commands are

<pre>module list                 # Show loaded modules
module unload <module>      # Unload a module
module keyword <keyword>    # All modules with a keyword (e.g. bam, alignment)</pre>

<div class="exercises">
Exercise<br/>
<ul>
<li>Load modules for bowtie2 and samtools.</li>
<li>Find and load modules for htseq-count,  macs and homer</li>
</div>

#### Alignment

For this workshop we're going to go straight to aligning the reads to the appropriate genome using bowtie and output an (indexed and sorted) bam file.   This is an alignment file and we will visualize this in the IGV browser and use it as input for the peak calling software. To align to a genome we need the fasta genome file specially indexed.   For the common genomes we have prepared indices on the regal filesystem. 

The indexed mouse mm10 genome files are here :

    /n/regal/informatics_public/ref/ucsc/Mus_musculus/mm10/chromFa.fa

#### The bowtie2 alignment program

We're going to use the bowtie2 alignment program.   For paired end reads we need to run it once on each pair of files (so 4 commands - one for each sample).  The command line syntax (including the module commands) is

    :::bash
    module load centos6/bowtie2-2.2.1
    module load centos6/samtools-0.1.19
    
    bowtie2  \
    -x <genome index>  \
    -1 <read1 fastqfile> \
    -2 <read2 fastqfile> \
    -p <cores> \
    > out.sam

This produces a file format called sam.   This has the nice property of being human readable but the downside of being rather large. We can use samtools to process the sam file into a smaller, compressed and indexed format called bam.    We also need to sort the alignment so it can be property indexed. We could run this as separate commands on the bowtie2 output as follows :

    :::bash
    samtools view -b -S out.sam > tmp.bam
    samtools sort tmp.bam out
    samtools index out.bam


This has the disadvantage of creating a lot of intermediate files that we won't use and take up space on the filesystem.  Instead we will pipe the output into the different commands using |

    :::bash
    bowtie2  \
    -x <genome index>  \
    -1 <read1 fastqfile> \
    -2 <read2 fastqfile> \
    -p <cores>  | samtools view -b -S - | samtools sort - out
    
    samtools index out.bam

Much neater yes?  This does have the disadvantage that if anything goes wrong in the samtools view or samtools sort part we have to do the whole alignment again but the tools are robust enough that that is rare. 

<div class="exercises">
Exercise<br/>

Construct a bowtie2 command line for one of the fastq files.   Do a test run (ctrl-C will interrupt a command) to check it starts running.
</div>

#### Constructing a Slurm submit script

Sadly things get a little more complicated still.    We can't just type in a command and run it on the login nodes (apart from testing).   We need to put our commands in something called an sbatch submit script. This script has a number of lines at the top that tell the scheduling system Slurm how and where to run the job.  Each time you use this script to submit a job you need to check the lines starting #SBATCH

    :::bash
    #!/bin/bash
    
    # Using $1 $2 etc in the command means you can call this script with an input file 
    # i.e. to submit you would do
    # sbatch [thisfile.sh](http://thisfile.sh) file1.dat file2.dat

    #SBATCH -J <jobname>
    #SBATCH -N 1                                 # Ensure that all cores are on one machine
    #SBATCH -n <cores>                           # Use <cores> cores for one job
    #SBATCH -t D-HH:MM                           # Runtime in D-HH:MM
    #SBATCH -p <partition>                       # Partition to submit to
    #SBATCH --mem=<mem>                          # Memory pool in Mb for all cores 
    #SBATCH -o <outfile>.%A.out                  # File to which STDOUT will be written
    #SBATCH -e <outfile>.%A.err                  # File to which STDERR will be written
    #SBATCH --mail-type=ALL                      # Type of email notification- BEGIN,END,FAIL,ALL
    #SBATCH --mail-user=michele.clamp@gmail.com  # Email to which notifications will be sent

    # Now your module commands - you can have multiple
    module load centos6/mymodule

    # Now your command (using input parameters $1 $2 etc)
    mycommand $1 $2 > $1.out

Note:  The $1 and $2 refer to arguments you put on the command line.      If this script was saved into a file myfile.sh we would submit this using the sbatch command as follows

    sbatch myfile.sh  file1.fast1 file2.fastq

i.e. $1 refers to the fist argument file1.fastq and $2 refers to the 2nd argument file2.fastq. This way we can have a submit script that can be reused on different files.

#### Summary of SLURM commands

The table below shows a summary of SLURM commands, along with an example. These commands are described in more detail below along with links to the SLURM doc site.

<table>
<tbody>
<tr>
<th>SLURM</th>
<th>SLURM EXAMPLE</th>
</tr>
<tr>
<td>Submit a batch serial job</td>
<td>sbatch</td>
<td><code>sbatch [runscript.sh]</code></td>
</tr>
<tr>
<td>Run a script interatively</td>
<td>srun</td>
<td><code>srun --pty -p interact -t 10 --mem 1000 /bin/hostname</code></td>
</tr>
<tr>
<td>Kill a job</td>
<td>scancel</td>
<td><code>scancel 999999</code></td>
</tr>
<tr>
<td>View status of queues</td>
<td>squeue</td>
<td><code>squeue -u akitzmiller</code></td>
</tr>
<tr>
<td>Check current job by id</td>
<td>sacct</td>
<td><code>sacct -j 999999</code></td>
</tr>
</tbody>
</table>

Also see [here](https://rc.fas.harvard.edu/docs/running-jobs.html) for the Research Computing page which has more info on the sbatch commands.

#### Constructing a bowtie2 submit script.

We are now going to combine out bowtie2/samtools command with the slurm script to produce a script that can be submitted to the cluster.  We've produced a template for this for you to modify.

    #!/bin/bash 

    # Call this script with 2 fastq file2 for reads 1 and 2.
    #
    # e.g. sbatch [thisscript.sh](http://thisscript.sh) myfile.R1.fastq  myfile.R2.fastq
    #
    # myfile.R1.fastq is referenced by the variable $1
    # myfile.R2.fastq is referenced by the variable $2 

    #SBATCH -J ATAC_Bowtie2 
    #SBATCH -N 1                                 # Ensure that all cores are on one machine
    #SBATCH -n <cores>                           # Use n cores for one job 
    #SBATCH -t 0-12:00                           # Runtime in D-HH:MM 
    #SBATCH -p general                           # Partition to submit to 
    #SBATCH --mem=<mb_ram>                       # Memory pool for all cores 
    #SBATCH -o <outfile>.%A.out                  # File to which STDOUT will be written 
    #SBATCH -e <outfile>.%A.err                  # File to which STDERR will be written 
    #SBATCH --mail-type=ALL                      # Type of email notification- BEGIN,END,FAIL,ALL 
    #SBATCH [--mail-user=<](mailto:--mail-user=michele.clamp@gmail.com)myemail>                # Email to which notifications will be sent 

    module load centos6/bowtie2-2.2.1 
    module load centos6/samtools-0.1.19 

    bowtie2 -x <indexfile> \ 
    -1 $1 \ 
    -2 $2 \ 
    -p <cores> | samtools view -b -S - |samtools sort - $1 

    samtools index $1.bam

<div class="exercises">
Exercise<br/>
Take this template script and modify it in nano so it has the following: 
<li>a suitable job name </li>
<li>runs on 32 cores (change this in 2 places) </li>
<li>runs in the general partition </li>
<li>uses 32Gb RAM </li>
<li>has a suitable output filename </li>
<li>sends email to yourself </li>
<li>the correct mm10 index file (see previous section)</li>
</div>

#### Submitting to slurm

An example slurm submission command is :

    :::bash
    sbatch  sbatch.bowtie2.sh Proximal-Femur-1_S1.R1.fastq Proximal-Femur-1_S1.R2.fastq

Test your submit script with this command and then construct and submit commands for the remaining 3 samples.

* use `squeue -u <username>` to check on your submissions 
* check the output and error files to see if it's running correctly. 

Note it's perfectly allowed to test your command on the login node just for a couple of seconds to check it does actually run before putting it in an sbatch script.

#### Alignment Files

There are pre-computed bam files in this directory : 

    /n/regal/informatics/workshops/Capellini/Output/ 
    
Make symlinks to your user directory (like you did for the fastq files) so you can work on them in your own directory.

    :::bash
    for i in ../../Output/*.bam ; do
       ln -s $i .

    for i in ../../Output/*.bam.bai; do
       ln -s $i .</pre>

### Alignment Statistics

Use samtools to get coverage stats for your alignment and then visualize the bam in IGV on your laptop.

    :::bash
    module load centos6/samtools-0.1.19
    
    for i in *.bam; do
      samtools flagstat $i

You can also use other modules such as bamtools or get gene level counts using htseq-count 

<div class="exercises">
Exercises<br/>
<ul>
<li>Run samtools on your bamfiles using a for loop </li>
<li> Find the bamtools module and load it </li>
<li> Run the bamtools command in a for loop on your bamfiles </li>
<li>(Extra:   Instead of using *.bam in your for loop use a find command instead) </li>
<li>(Extra:   Find and load the htseq-count module and run it on a bam file.</li>
</ul>
</div>

#### Visualize BAMS in IGV

Transfer the bam files  and the bam index files to your laptop (or to a directory that you can mount on your laptop).   The full bam files are very big so we've generated a smaller slice of the bam files for a single chromosome in

    /n/regal/informatics/workshops/Capellini/Output/chr12bam.tgz

From your laptop do :

    :::bash
    scp -r username@login.rc.fas.harvard.edu:/n/regal/informatics/workshops/Capellini/Output/chr12bam.tgz .

Note:  to split a set of bam files by chromosome (ie. by reference) do

    :::bash
    module load bio/bamtools-2.2.0
    for i in *.bam; do
      bamtools split -in $i -reference
    done

#### Install and run IGV

Registration and installation instructions are here https://www.broadinstitute.org/igv/ Make sure the bam and the bam.bai files are in the same directory and load the bam file.   Zoom in until you can see reads - do the tracks look how you expect them to?

### MACs Peak Calling

We’re going to use MACS for peak calling which is a commonly used program and often comes out top in a comparison of methods.  The publication describing the algorithm is [here](http://www.genomebiology.com/content/9/9/R137). Paramount in any peak calling method is the use of controls as reads are not distributed randomly across the genome and can generate many false positives. MACS empirically models the shift size of ChIP-Seq tags, and uses it to improve the spatial resolution of predicted binding sites. MACS also uses a dynamic Poisson distribution to effectively capture local biases in the genome.

#### MACS Command Line

An example of a macs command line with commonly used options is below

    :::bash
    macs2 callpeak  \
    --name macs1
    --treatment file1.bam file2.bam \
    --control   ctrl1.bam ctrl2.bam  \
    --outdir macs1 \
    --format BAM \
    --name macs1 \
    --pvalue 1e-5 \
    --gsize hs \
    --tsize 75

<table>
<tr><td nowrap><code>--treatment</code></td><td>The condition alignments in sam format.  Multiple files can be specified separated by spaces.</td></tr>
<tr><td nowrap><code>--control</code></td><td>The control alignments in sam format</td></tr>
<tr><td nowrap><code>--pvalue</code></td><td>Use 1e-5 at least</td></tr>
<tr><td nowrap><code>--name</code></td><td>The output directory where the results will go</td></tr>
<tr><td nowrap><code>--gsize</code></td><td>Genome size.  Use hs for human or the actual base pair size</td></tr>
<tr><td nowrap><code>--format</code></td><td>SAM (can also use BAMs and other formats)</td></tr>
<tr><td nowrap><code>--tsize</code></td><td>Read length</td></tr>
<tr><td nowrap><code>--bdg</code></td><td>Create a bed graph  (optional)</td></tr>
<tr><td nowrap><code>--wig</code></td><td>Create a wig file  (optional)</td></tr>
</table>

The only essential parameter is the `--treatment` one but I recommend you include the others (up to tsize) as a record of how you ran the program. More parameters and their descriptions are at [https://github.com/taoliu/MACS/blob/macs_v1/README.rst](https://github.com/taoliu/MACS/blob/macs_v1/README.rst)

#### Suggestions for MACS sbatch files

<table>
<tr><td>Running time</td><td>~1 hour per sample</td></tr>
<tr><td>Memory needs</td><td>~2G</td></tr>
<tr><td>Module</td><td><a href="https://portal.rc.fas.harvard.edu/apps/modules/macs2/2.1.0.20140616-fasrc01">macs2/2.1.0.20140616-fasrc01</a></td></tr>
<tr><td>Queue</td><td>serial_requeue usually schedules pretty quickly but jobs can get cancelled.</td></tr>
</table>

#### Output files

The main output file is  NAME_peaks.narrowpeaks which is a  BED format file which contains the peak locations.   This will be the file you use as input for the motif analysis.

You can load it to UCSC genome browser or Affymetrix IGB software. The 5th column in this file is the `-10*log10` pvalue of peak region.

Other output formats are described at [https://github.com/taoliu/MACS/blob/macs_v1/README.rst#output-files](https://github.com/taoliu/MACS/blob/macs_v1/README.rst#output-files)

#### Exercise:   Peak Calling

<div class="exercises">
<ol>
<li> Construct an sbatch submission script for your treatment and control files and submit to slurm. </li>
<li>Construct a reverse sbatch submission script using the previous treatment files as controls and the previous control files as treatment. </li>
<li>Use squeue -u <username> to check the jobs are submitted and either pending or running.   Look at the error files if if fails to run. </li>
<li>Look at the output   (see <code>/n/regal/informatics/workshops/Capellini/Output/macs2*</code>  for pre-computed output) </li>
<li> Download the peaks file  to your laptop </li>
<li> Load up your bam file and peaks file into IGV </li>
<li> Zoom and scroll through some peaks.   </li>
</ol>

Look at some high scoring peaks and some low scoring peaks.  Is this what you expect?   Is there anything different (bar the obvious) between the high scoring and the low scoring peaks.
</div>

## Peak Annotation Using Homer

[http://homer.salk.edu/homer/introduction/basics.html](http://homer.salk.edu/homer/introduction/basics.html)

One of the many things that Homer can do is take a set of regions in a genome (our peaks for example) and annotate them with respect to where they lie compared to genes.  Each peak is assigned to a gene using the closest TSS (transcription start site).

The command looks like

    annotatePeaks.pl \
    <peakfile> \
    <genome>   \
    > <outputfile>


You can load the output file into excel if you want to look at your individual peaks later.

[Output columns](http://homer.salk.edu/homer/ngs/annotation.html) .

    Description of Columns:

    1.  Peak ID
    2.  Chromosome
    3.  Peak start position
    4.  Peak end position
    5.  Strand
    6.  Peak Score
    7.  FDR/Peak Focus Ratio/Region Size
    8.  Annotation (i.e. Exon, Intron, ...)
    9.  Detailed Annotation (Exon, Intron etc. + CpG Islands, repeats, etc.)
    10.  Distance to nearest RefSeq TSS
    11.  Nearest TSS: Native ID of annotation file
    12.  Nearest TSS: Entrez Gene ID
    13.  Nearest TSS: Unigene ID
    14.  Nearest TSS: RefSeq ID
    15.  Nearest TSS: Ensembl ID
    16.  Nearest TSS: Gene Symbol
    17.  Nearest TSS: Gene Aliases
    18.  Nearest TSS: Gene description
    19.  Additional columns depend on options selected when running the program.


<div class="exercises">
Exercise<br/>
<ul>
<li>Construct a command line for the annotatePeaks.pl command using one of the peaks files and the mm10 genome. </li>
<li>Run this (you don't need a slurm script for this) and generate an output file </li>
<li>Either look at the file on the cluster or copy to your laptop and view in excel.   Is this output useful? Example output is in :  /n/regal/informatics/workshops/Capellini/Output/macs2/macs2_all_cmdline_peaks.annotatePeaks and /n/regal/informatics/workshops/Capellini/Output/macs2/_rev/macs2_all_cmdline.rev_peaks.annotatePeaks </li>
<li>Compare the two files roughly (number of peaks, max/min score of peaks, chromosomal distribution of peaks).   Do you expect the differences?</li>
</ul>
</div>

#### GO Term Enrichment Using Homer

So right now we have a list of regions and the genes they are associated with.   We could just take the top few high scoring genes as a candidate list for further exploration.    However it's rare that the effect we're investigating is down to one or two highly influential genes.   Usually a group of genes is responsible.    To investigate this we can look at the list of genes and see whether their assigned function is statistically significant. This can be run for mouse using the findGO.pl command in the homer suite.  The command line looks like

    :::bash
    findGO.pl <file containing entrez ids> <genome> <outdir>

The entrez ids are conveniently found in the annotatePeaks output in column 12 and you can extract them using the following awk command  

    :::bash
    awk -F'\t' 'NR > 1 {print $12}'  <annotatepeaks_outfile> > <entrezidsfile>

This can be run on the command line - no need for sbatch.

Using this file you can then run the GO enrichment command

    :::bash
    findGO.pl <entrezidsfile> mouse <outdir>

The output gets written to a directory which contains (amongst other things) an HTML file.  I generally download the whole directory to my laptop/desktop using scp and open the HTML file in a browser. 

<div class="exercises">
Exercise<br/>
<ul>
<li>Using your (or the supplied) annotatePeaks output file (either one) extract the entrezids into a file</li>
<li>Run findGO.pl on this file using the mouse genome and write to a directory findGO </li>
<li>Copy the findGO directory to your laptop and open the html file, e.g. I would do
<div class="highlight"><pre>scp -r mclamp@sandy2:/n/regal/informatics/workshops/Capellini/Output/macs2/findGO findGO
open findGO/geneOntology.html</pre></div></li>
</ul>
Look at the enriched GO terms - do they make sense? 

<ul>
<li>Do this for the other (reversed) set of annotated peaks and look at the enriched GO terms.</li>
<li>Does it make sense to combine the two sets of data and analyse together? </li>
<li>Do we want to filter our peaks first?</li>
</ul>
</div>

## Extras:

### Motif Finding using Homer

[http://homer.salk.edu/homer/introduction/basics.html](http://homer.salk.edu/homer/introduction/basics.html) Motif finding has many pitfalls.

*   The motifs are small and common in the genome,
*   The motifs are not always perfectly described in the database
*   There is a lot of overlap beween motifs (degeneracy)
*   Methods look for differential enrichment between a condition and control set of sequences.
*   Controls are hard to get right and need thought.

By default HOMER will use confident, non-regulated promoters as background when analyzing promoters, and sequences in the vicinity of genes for ChIP-Seq analysis (i.e. from –50kb to +50kb).  In each case sequences are matched for their GC content to avoid bias from CpG Islands.

Homer motif finding usage
[http://homer.salk.edu/homer/ngs/peakMotifs.html](http://homer.salk.edu/homer/ngs/peakMotifs.html)

    :::bash
    findMotifsGenome.pl <mypeakfile>  \
    <mygeneome|mygenomefile>  \
    <outputdir> \
    -size <regionsize> \
    -p <threads> \
    -len <len1>,<len2>,<len3>
    
<table>
<tr><td><code>mypeakfile</code></td><td>the peak file from the macs output (.narrowpeak)</td></tr>
<tr><td><code>mygenome</code></td><td>hg19,mm9</td></tr>
<tr><td><code>mygenomefile</code></td><td>hg19.fa mm9.fa (it will preprocess these)</td></tr>
<tr><td><code>outputdir</code></td><td>output directory</td></tr>
<tr><td><code>-size</code></td><td>the size of the region (default is the whole peak)  (start with 50 for finding the main factor and work out)</td></tr>
<tr><td><code>-p</code></td><td>number of threads</td></tr>
<tr><td><code>-len</code></td><td>the lengths of the motifs to consider separated by commas</td></tr>
<tr><td><code>-S</code></td><td>number of motifs to find.  The default is 25 which the docs say is pretty high already</td></tr>
<tr><td><code>-dumpFasta</code></td><td>dumps the fasta of the regions</td></tr>
</table>

### Module commands needed

    :::bash
    module load homer

### Hands on

1. Construct and submit an sbatch file to find motifs on the MACS peaks for a motif of length 7 in a region of 50bp. Use multiple threads (16 is good),  1Gb memory and 30 mins run time.

2. Construct and submit an sbatch file to find motifs on the MACS peaks for motifs of lengths 7 in a region of 200bp.  Use multiple threads, 1Gb memory and 30 mins run time.<br/>
Either wait for the jobs to run or look at the pre-run output in

    `/n/regal/informatics/workshops/ChIP-Seq/Output/macs2_all/macs2_all_peaks.narrowPeak_homer_50bp_7mer` <br/><br/>
or for the 200bp region

    `/n/regal/informatics/workshops/ChIP-Seq/Output/macs2_all/macs2_all_peaks.narrowPeak_homer_100bp_7mer`
    
3. Transfer the output directory to your laptop and look at the html output.   Considering the protein used in this experiment are the results understandable?   Do the differences between the 50bp and 200bp region make sense.

4. How would you construct an sbatch file to search for motifs of multiple lengths?


