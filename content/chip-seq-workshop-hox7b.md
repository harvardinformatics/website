Title: Chip-Seq Workshop
Date: 2016-05-09
Category: Tutorials
Author: Michele Clamp
Tags: Chip-Seq, Workshop
Summary: This tutorial provides a workflow for Chip-Seq alignment, peak finding and motif analysis

## Notes

The input data files are available at:

    /n/regal/informatics/workshops/ChIP-Seq/Data
    
The sbatch files shown below are available at:

    /n/regal/informatics/workshops/ChIP-Seq/Slurm

The output files for the intermediate steps are available at:

    /n/regal/informatics/workshops/ChIP-Seq/Output
  
## Prerequisites

* Ability to login to odyssey (ssh, openauth)
* Knowledge of basic command line commands (ls, cd, cat, less, nano, mkdir, |, >)
* Knowledge of basic script loops
* Knowledge of slurm batch scripts
* Ability to transfer files to/from the cluster
* Have IGV installed on your laptop. (https://www.broadinstitute.org/software/igv/download)

## Chip-Seq Overview

Chip-seq, or Chromatin Immunoprecipitation Sequencing, is used to probe DNA-protein interactions. First, DNA in whole cells is reversibly cross-linked to bound proteins with a cross-linking agent such as formaldehyde. Next, the DNA is isolated and sonicated to produce small fragments. Fragments containing the protein of interest are isolated via the addition of a protein-specific antibody and chromatin immunoprecipitation. Crosslinks are reversed and the isolated DNA is purified. Finally, strand-specific adapters are ligated and the library is amplified.

<figure>
	<a class="img" href="/images/chipseq/overview.png">
    		<img class="img-responsive" src="/images/chipseq/overview.png"></img>
	</a>
    <figcaption>Figure 1. ChIP-Seq Overview [Park 2009] </figcaption>
</figure>


## Experimental Design

When starting to think about a ChIP-Seq sequencing experiment several things need to be considered.  These are :

*  number of replicates
*  suitable controls
*  number and type of reads

Guidelines on these questions are addressed in the following publications :

* [ENCODE guidelines](http://www.ncbi.nlm.nih.gov/pubmed/22955991) Landt et al 2010
* [Practical Guidelines for the Comprehensive Analysis of ChIP-Seq data](http://www.ncbi.nlm.nih.gov/pmc/articles/pmid/24244136/) Bailey et al 2013

### Sequencing depth
The number of reads needed depends on the type and frequency of peak expected (sharp/broad,  common/rare,   long motif/short motif) and also the size of the genome.

Rough guidelines are as follows:

#### Mammalian cells
Sharp peaks (TFs)

* 10 million uniquely mapped (with 80% distinct reads) per replicate


Broad peaks (chromatin mods e.g. H3H27me3 M3K36me3)

* 12-20 million uniquely mapped (with 80% distinct reads) per replicate


#### Flies/Worms
Sharp peaks (TFs)

*  2 million uniquely mapped (with 80% distinct reads) per replicate

Broad peaks (chromatin mods e.g. H3H27me3 M3K36me3)

*  5-10 million uniquely mapped (with 80% distinct reads) per replicate

### Controls

These are essential to reduce false positives as there is always some signal on open chromatin and there is a correlation between the number of reads found in protein-bound regions and non protein bound regions.  These artefacts originate from a number of sources including :

 * Copy number variation
 * Non-uniform fragmentation
 * Non-specific pull-down
 * Incorrect mapping of repetitive genomic regions
 * GC sequencing bias (http://beads.sourceforge.net [Cheung et al 2011])
 * Larger genomes have more problems.
 
There are 2 different types of controls used for ChIP-Seq

 * Input DNA (sonicated chromatin - no IP)
 * DNA with a different antibody (often IgG)
 
 Some people prefer input DNA as IgG yields can be much lower.  This affects the statistics when we come to detecting significant peaks in the data (see the MACS paper Zhang and Liu 2008).
 
### Replicates

The good news here is that, unlike RNA-Seq, more than 2 replicates does not significantly increase the number of targets.

### Length and type of reads

Single ended 50bp reads are good for ChIP-Seq.    We don't need the extra cost of paired-end reads and 50bp is long enough for mapping to a reference genome.


## Outline of ChIP-Seq Sequencing and Analysis Process

<figure>
	<a class="img" href="/images/chipseq/analysis_workflow.png">
    		<center><img class="img-responsive" src="/images/chipseq/analysis_workflow.png"></img></center>
	</a>
    <figcaption>Figure 3.  ChIP-Seq Analysis Workflow</figcaption>
</figure>


Strand-specific reads help to distinguish signal from noise, as true signal will consist of two peaks of equal magnitude, but on opposite strands. This is due to the fact that reads can originate from either end of DNA fragments, and DNA fragments are usually much longer than reads. This characteristic peak profile is shown in Figure 3.

<figure>
	<a class="img" href="/images/chipseq/read_peak_profile.png">
    		<center><img class="img-responsive" src="/images/chipseq/read_peak_profile.png"></img></center>
	</a>
    <figcaption>Figure 2. ChIP-Seq Read Peak Profile [Park 2009]</figcaption>
</figure>


## Workshop dataset
[http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47164](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47164)

Examination of chromatin binding sites of HOXB7 in BT-474 breast cancer cell line using ChIP-seq. Four parallel IgG samples were sequenced, merged together and used as a control data set. Two parallel HOXB7 ChIP samples were sequenced and merged for each replicate, AF1 and AF2. Both HOXB7 ChIP replicates (AF1 and AF2) contained approximately the same amount of reads as the merged IgG control data set.

This is a recent dataset with replicates and controls for a single transcription factor which should give a healthy number of narrow well-defined peaks.  Ideal for a workshop!


The associated publication is <a href="http://onlinelibrary.wiley.com/doi/10.1002/ijc.29616/abstract">http://onlinelibrary.wiley.com/doi/10.1002/ijc.29616/abstrac</a>t (<a href="http://onlinelibrary.wiley.com/doi/10.1002/ijc.29616/epdf">pdf</a>)

Treatment and control sequence files are :

<table cellpadding="3">
<tbody>
<tr>
<td valign="top">
<div><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1145931">GSM1145931</a></div></td>
<td valign="top">
<div>HOXB7_AF1a   SRR866428</div></td>
</tr>
<tr>
<td valign="top">
<div><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1145932">GSM1145932</a></div></td>
<td valign="top">
<div>HOXB7_AF1b    SRR866429</div></td>
</tr>
<tr>
<td valign="top">
<div><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1145933">GSM1145933</a></div></td>
<td valign="top">
<div>HOXB7_AF2a    SRR866430</div></td>
</tr>
<tr>
<td valign="top"><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1145934">GSM1145934</a></td>
<td valign="top">
<div>HOXB7_AF2b   SRR866431</div>
<div></div></td>
</tr>
<tr>
<td valign="top">
<div><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1145935">GSM1145935</a></div></td>
<td valign="top">
<div>IgGa     SRR866432</div></td>
</tr>
<tr>
<td valign="top"><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1145936">GSM1145936</a></td>
<td valign="top">IgGb      SRR866433</td>
</tr>
<tr>
<td valign="top"><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1145937">GSM1145937</a></td>
<td valign="top">
<div>IgGc      SRR866434</div></td>
</tr>
<tr>
<td valign="top">
<div><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1145938">GSM1145938</a></div></td>
<td valign="top">
<div>IgGd      SRR866435</div>
<div></div></td>
</tr>
</tbody>
</table>

The raw fastq files are available on the cluster in

    /n/informatics/workshops/ChIP-Seq/Data

## Cluster login and setup

1. If you haven't already done so log into the cluster and cd into the workshop directory
<div class="codehilite"><pre>
ssh login.rc.fas.harvard.edu
cd /n/regal/informatics/workshops/ChIP-Seq/Users
</pre></div>

2. Make your own user directory if you haven't already and cd into it
<div class="codehilite"><pre>
mkdir myusername
</pre>
</div>

3. Make a set of directories where you will store your data, script files and output.   Analyses like these can end up generating a lot of files so keeping things tidy pays off.

<div class="codehilite"><pre>
cd myusername
mkdir Data
mkdir Output
mkdir Slurm
</pre>
</div>

4. Now we're going to change into the Data directory and do a fancy bash loop to link to the datafiles. (We could copy them over but this saves disk space).

<div class="codehilite"><pre>
cd Data
for i in /n/regal/informatics/workshops/ChIP-Seq/Data/*.fastq ; do
   ln -s $i .
done

ls -l 
cd ..
</pre>
</div>

## Pre-processing and Alignment

### Adapter and quality trimming

These steps are good practice but it could be argued that they aren't essential for this type of analysis.  Unlike assembly or variant calling we're not interested in the exact alignment of bases but are only interested in the read placement on the genome.  If the alignment is slightly off by a few bases it's not going to ruin things.

We're going to submit a job to the cluster using the scheduling software Slurm.   Let's take a look at a template Slurm submit script (/n/regal/informatics/workshops/ChIP-Seq/Slurm/sbatch.template.sh)

<div class="codehilite">
<pre>
#/bin/bash
# Using $1 in the command means you can call this script with an input file
#
# Change the &lt;nnnn&gt; strings to something sensible before submitting
#
#SBATCH -J &lt;jobname&gt;
#SBATCH -N 1                   # Ensure that all cores are on one machine
#SBATCH -n &lt;n&gt;                 # Use n cores for the job
#SBATCH -t &lt;n-nn:nn&gt;           # Runtime in D-HH:MM
#SBATCH -p &lt;queuename&gt;         # Partition to submit to
#SBATCH --mem=&lt;n&gt;              # Memory pool for all cores in Mb (see also --mem-per-cpu)
#SBATCH -o &lt;outfile&gt;.%A.out    # File to which STDOUT will be written (%A is replaced by the jobid)
#SBATCH -e &lt;outfile&gt;.%A.err    # File to which STDERR will be written (%A is replaced by the jobid)
#SBATCH --mail-type=&lt;type&gt;     # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=&lt;myemail@harvard.edu&gt; # Email to which notifications will be sent

source new-modules.sh
module load <mymodule>

# The command can use the input parameter $1 here.  Extra command line arguments can be used with $2 $3 etc

some_command_here  $1 > $1.out
</pre>
</div>

## Summary of SLURM commands
The table below shows a summary of SLURM commands, along with LSF equivalents and an example. These commands are described in more detail below along with links to the SLURM doc site.

<table>
<thead>
<tr>
<th></th>
<th>SLURM</th>
<th>SLURM Example</th>
</tr>
</thead>
<tbody>
<tr>
<td>Submit a batch serial job</td>
<td>sbatch</td>
<td><code>sbatch runscript.sh</code></td>
</tr>
<tr>
<td>Run a script interactively</td>
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

Also see <a href="https://rc.fas.harvard.edu/docs/running-jobs.html">here</a> for the Research Computing page which has more info on the sbatch commands.

There are example Slurm scripts in 

<pre>/n/regal/informatics/workshops/ChIP-Seq/Slurm</pre>

### Exercises 1

* Make a trimmomatic directory in your Output directory and change into it.  This is good practice as it keeps things neat and tidy.
* Copy over the example trimmomatic slurm script to your Slurm directory
* Edit it and change the email address.  Check the parameters look good.
* Test the script on a fastq file using 

<pre> bash sbatch.Trimmomatic_SE.sh myfile.fastq</pre>

  If this doesn't exit immediately you probably have the right syntax so you're good to submit to slurm.
* Submit a slurm job for each of the fastq files using a for loop.  The submit syntax is: 

<pre>sbatch sbatch.Trimmomatic_SE.sh myfile.fastq</pre>

* Use the squeue command to check on your jobs.

When you have your jobs successfully in the queue we'll let them stay there.   They take too long to run for us to wait for them so we'll take a look at the pre-computed output.

We've provided output files for each of the steps in 

  <pre>/n/regal/informatics/workshops/ChIP-Seq/Output</pre>

For the following steps use these as input files where needed.

### Alignment

We're going to use bowtie2 to align our reads to the genome.  As well as the input read files bowtie needs an indexed genome which can be found here:

<pre>/n/regal/informatics_public/ref/ucsc/Homo_Sapiens/hg19/chromFa</pre>

An example slurm submission script is here 

<pre>/n/regal/informatics/workshops/ChIP-Seq/Slurm/sbatch.bowtie2_SE.sh</pre>

<pre>
#!/bin/bash

# Call this script with an input refrence genome and a single fastq file.
#
#  e.g. sbatch sbatch.bowtie2_SE.sh /n/regal/informatics_public/ref/ucsc/Mus_musculus/mm10/chromFa myfile.fq
#
# The aligned output will go to myfile.bam
#
# Change the -o line if you need to change the output file

#SBATCH -J Bowtie2_SE
#SBATCH -N 1                                 # Ensure that all cores are on one machine
#SBATCH -n 16                                # Use 16 cores for one job
#SBATCH -t 0-4:00                            # Runtime in D-HH:MM
#SBATCH -p serial_requeue                    # Partition to submit to
#SBATCH --mem=32000                          # Memory pool for all cores
#SBATCH -o bowtie2.%A.out                    # File to which STDOUT will be written
#SBATCH -e bowtie2.%A.err                    # File to which STDERR will be written
#SBATCH --mail-type=ALL                      # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com  # Email to which notifications will be sent

module purge

source new-modules.sh

module load bowtie2
module load samtools

module list

# This gets the stub of the filename without .gz .fq or .fastq

j=$( basename $2 )
j=${j%.gz}
j=${j%.fastq}
j=${j%.fq}

# Now run the command

bowtie2 -x $1 \
-U $2 \
-p 16 \
 > $j.sam

samtools view -b -S $j.sam | samtools sort - $j   # Convert to bam and sort
samtools index $j.bam                             # Index the bam
</pre>

This will produce a .bam alignment file and a .bam.bai index file.    The job error file will contain some statistics on the alignment.


### Exercises 2

 * Make a bowtie2 directory in your Output directory (tidiness first!) and cd into it.
 * Copy over the bowtie2_SE slurm script from the workshop Slurm directory and edit it to change the email.
 * Test the script using 

<pre>bash sbatch.bowtie2_SE.sh myfile.fastq</pre>

You can use the raw *.fastq files in the Data directory or the trimmed ones in the workshop output directory

 * Once your script is running submit the job to slurm for each of the fastq files. 

<pre>sbatch sbatch.bowtie2_SE.sh myfile.fastq</pre>

For extra points use a for loop to submit all the jobs.

 * Use squeue to check on the progress of your jobs.

## MACS Peak Calling

We’re going to use MACS2 for peak calling which is a commonly used program and often comes out top in a comparison of methods.  The publication describing the algorithm is [here](http://www.genomebiology.com/content/9/9/R137).

Paramount in any peak calling method is the use of controls as reads are not distributed randomly across the genome and can generate many false positives.

In addition MACS empirically models the shift size of ChIP-Seq tags, and uses it to improve the spatial resolution of predicted binding sites. MACS also uses a dynamic Poisson distribution to effectively capture local biases in the genome.


<figure>
	<a class="img" href="/images/chipseq/macs_poisson.png">
    		<center><img class="img-responsive" src="/images/chipseq/macs_poisson.png"></img></center>
	</a>
    <figcaption>Figure 4. Using dynamic Poisson model for peak calling</figcaption>
</figure>


## MACS Fragment Size Estimation

We sequence only one end of the sonicated and protein bound fragment.   The actual protein binding site will be somewhere along the fragment length.  To more precisely locate the binding site we need to shift the alignment position by half the mean fragment length.

MACS randomly samples 1000 high quality (high fold enrichment) peaks, separates their forward and reverse tags and aligns them by their mid point.   The distance between the modes is labeled the fragment size d.

<figure>
	<a class="img" href="/images/chipseq/macs_fragment_estimation.png">
    		<center><img class="img-responsive" src="/images/chipseq/macs_fragment_estimation.png"></img></center>
	</a>
    <figcaption>Figure 5.  Using read alignments to estimate fragment size</figcaption>
</figure>


## MACS2 Command Line
An example of a macs command line with commonly used options is below
<pre>
    macs2 callpeak  \
    --name macs_run1
    --treatment file1.bam file2.bam \
    --control ctrl1.bam ctrl2.bam \    # Optional (but desirable)
    --outdir macs_run1 \
    --format BAM \                     # Optional
    --pvalue 1e-5 \                    # P value
    --gsize hs \                       # Genome size
</pre>

We can also run MACS2 with a qvalue (FDR) threshold as follows

<pre>
    macs2 callpeak  \
    --name macs_run1
    --treatment file1.bam file2.bam \
    --control ctrl1.bam ctrl2.bam \    # Optional (but desirable)
    --outdir macs_run1 \
    --format BAM \                     # Optional
    --qvalue .01 \                     # FDR
    --gsize hs \                       # Genome size
</pre>


<table>
<tr><td nowrap><code>--treatment</code></td><td>The condition alignments in sam format.  Multiple files can be specified separated by spaces.</td></tr>
<tr><td nowrap><code>--control</code></td><td>The control alignments in sam format</td></tr>
<tr><td nowrap><code>--pvalue</code></td><td>Use 1e-5 at least</td></tr>
<tr><td nowrap><code>--qvalue</code></td><td>FDR adjusted pvalue</td></tr>
<tr><td nowrap><code>--name</code></td><td>The output directory where the results will go</td></tr>
<tr><td nowrap><code>--gsize</code></td><td>Genome size.  Use hs for human or the actual base pair size</td></tr>
<tr><td nowrap><code>--format</code></td><td>SAM (can also use BAMs and other formats)</td></tr>
<tr><td nowrap><code>--tsize</code></td><td>Read length</td></tr>
<tr><td nowrap><code>--bdg</code></td><td>Create a bed graph  (optional)</td></tr>
<tr><td nowrap><code>--wig</code></td><td>Create a wig file  (optional)</td></tr>
</table>

The only essential parameter is the `--treatment` one but I recommend you include the others (up to tsize) as a record of how you ran the program.

More parameters and their descriptions are at

[https://github.com/taoliu/MACS/blob/master/README.rst](https://github.com/taoliu/MACS/blob/master/README.rst)

## MACS2 Slurm Submission

### Suggestions for MACS sbatch files
<table><tbody>
<tr><td>Running time </td><td>~1 hour per 10 million aligned reads</td></tr>
<tr><td>Memory needs</td><td>~2G per 10 million reads</td></tr>
<tr><td>Module</td><td><a href="https://portal.rc.fas.harvard.edu/apps/modules/macs2/2.1.0.20140616-fasrc01">macs2/2.1.0.20140616-fasrc01</a></td></tr>
<tr><td>Queue</td><td>serial_requeue usually schedules pretty quickly but jobs can get cancelled.</td></tr>
</tbody></table>

## Output files

The main output file is  NAME_peaks.narrowpeaks which is a  BED format file which contains the peak locations.   This will be the file you use as input for the motif analysis.

You can load it to UCSC genome browser or IGV. The 5th column in this file is the -10*log10pvalue of peak region.

Other output formats are described at
[https://github.com/taoliu/MACS/blob/macs_v1/README.rst#output-files](https://github.com/taoliu/MACS/blob/macs_v1/README.rst#output-files)

## Exercises 3 - Peak Calling and visualization

1. Make a macs2 directory in your user space for the macs2 output.
2. Copy over and edit the macs2 Slurm script from the workshop directory.
3. Construct an sbatch submission script for the 4  treatment and  4 control files.  Test on the command line and then submit. 
4.  Look at the output   (see /n/regal/informatics/workshops/ChIP-Seq/Output/macs2/ for pre-computed output)
    <ul><li>Download the peaks file, wig or bdg files and your bam files and indices to your laptop</li>
	   <li>Load up your bam files/wigs etc into IGV</li>
	   <li>Zoom and scroll through some peaks.   Look at some high scoring peaks and some low scoring peaks.  Is this what you expect?   Is there anything different (bar the obvious) between the high scoring and the low scoring peaks.</li><li>Extra:  The output files AF1.bed, AF2.bed from the HOXB7 paper are in the output directory.  Load these into IGV  and compare to your peaks.   How similar are they?   Can you write a command line to find the intersection between your results and the paper's results?</li><li>Extra:  Using the bedtools module  `bedtools2/2.25.0-fasrc01` use the bedtools intersect command to compare the AF1.bed and AF2.bed files to your results.
	   </li>
	</ul>
    

Hint:   To read the help info for a bedtools command do the following:

<pre>  bedtools intersect 2>&1 | more</pre>

Bedtools prints the help text to stderr (there are two output streams stdout and stderr).  By default `|more` will only redirect the stdout.   Using the `2>&1` will direct the stderr (2) output into the stdout stream (1)

## Motif Finding using Homer

[http://homer.salk.edu/homer/introduction/basics.html](http://homer.salk.edu/homer/introduction/basics.html)

Motif finding has many pitfalls.

*  The motifs are small and common in the genome,
*  The motifs are not always perfectly described in the database
*  There is a lot of overlap beween motifs (degeneracy)
*  Methods look for differential enrichment between a condition and control set of sequences.
*  Controls are hard to get right and need thought.

By default HOMER will use confident, non-regulated promoters as background when analyzing promoters, and sequences in the vicinity of genes for ChIP-Seq analysis (i.e. from –50kb to +50kb).  In each case sequences are matched for their GC content to avoid bias from CpG Islands.

Homer motif finding usage
[http://homer.salk.edu/homer/ngs/peakMotifs.html](http://homer.salk.edu/homer/ngs/peakMotifs.html)

<pre>
	findMotifsGenome.pl &lt;mypeakfile&gt; \
	    &lt;mygeneome|mygenomefile&gt;  \
	    &lt;outputdir&gt; \
	    -size &lt;regionsize&gt; \
	    -p &lt;threads&gt; \
	    -len &lt;len1&gt;,&lt;len2&gt;,&lt;len3
</pre>
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

<div class="codehilite"><pre>module load homer</pre></div>

### Exercises 4 

1. Construct and submit an sbatch file to find motifs on the MACS peaks for a motif of length 7 in a region of 50bp. Use multiple threads (16 is good),  1Gb memory and 30 mins run time.

2. Construct and submit an sbatch file to find motifs on the MACS peaks for motifs of lengths 7 in a region of 200bp.  Use multiple threads, 1Gb memory and 30 mins run time.

Either wait for the jobs to run or look at the pre-run output in

<div class="codehilite"><pre>/n/regal/informatics/workshops/ChIP-Seq/Output/macs2/macs2_multi/macs2_multi_peaks.narrowPeak_homer_50bp_7mer</pre></div>
	
or for the 200bp region

<div class="codehilite"><pre>/n/regal/informatics/workshops/ChIP-Seq/Output/macs2_multi/macs2_multi_peaks.narrowPeak_homer_100bp_7mer</pre></div>

3. Transfer the output directory to your laptop and look at the html output.   Considering the protein used in this experiment are the results understandable?   Do the differences between the 50bp and 200bp region make sense.

4. How would you construct an sbatch file to search for motifs of multiple lengths?

## Using Homer to search for a particular motif

The motif reported in the paper is GCCNGGC.  How can we search for this motif in out data?


1. First we generate a motif file for our motif (no need for a slurm script for this)

<div class="codehilite"><pre>seq2profile.pl <consensus> [# mismatches] [name] &gt; output.motif</pre></div>

2. Then we use that motif file to search our peaks file

<div class="codehilite"><pre>findMotifsGenome.pl  <peakfile>  \ 
    	&lt;genome&gt; \ 
        &lt;outputdir&gt; \
        -size &lt;size&gt; \
        -p &lt;threads&gt; \
       -find &lt;motiffile&gt;
</pre></div>

### Exercises 5

We need to generate a motif file but we don't have a consensus.  Let's run with a fake consensus to generate a template.


<div class="codehilite"><pre>seq2profile.pl GCCAGGC  1 HOXB7 > HOXB7.motif</pre></div>

Now we can edit the HOXB7.motif file to edit the probabilites to make the middle character 0.25,0.25,0.25,0.25

(Edit the A position to reflect the correct probabilities)
<div class="codehilite"><pre>
    >GCCAGGC    HOXB7    2.76827819473531
    0.001    0.001    0.997    0.001
    0.001    0.997    0.001    0.001
    0.001    0.997    0.001    0.001
    0.997    0.001    0.001    0.001
    0.001    0.001    0.997    0.001
    0.001    0.001    0.997    0.001
    0.001    0.997    0.001    0.001
</pre></div>

Then we can run findMotifsGenome.pl with this file to find candidates.


<div class="codehilite"><pre>findMotifsGenome.pl  macs2_peaks.narrowPeak hg19 macs2_homer_200_HOXB7 -size 200 -find HOXB7.motif  -p 16 &gt; macs2_homer_200_HOXB7.out</pre></div>


### Exercises 6

1.  Run the seq2profile.pl command on the login node (no need for a slurm command for this).

2.  Edit the resulting file to have the correct probability distribution for the middle base.

3. Construct a slurm batch file for the findMotifsGenome.pl and submit to the cluster.  Use a run time of ~1hour and 1Gb RAM.

To save time the output is in :

<div class="codehilite"><pre>/n/regal/informatics/workshops/ChIP-Seq/Output/macs2_all/macs2_homer_200_HOXB7</pre></div>

How do your results compare to the results in the paper?     What do you think is different?
