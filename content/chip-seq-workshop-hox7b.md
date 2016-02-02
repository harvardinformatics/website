Title: Chip-Seq Workshop - Hox7b
Date: 2015-01-01
Category: Tutorials
Author: Michele Clamp
Tags: Chip-Seq, Workshop
Summary: This tutorial provides useful tips for Chip-Seq experiment analysis.

## Notes
The sbatch files shown below are available at:

    /n/regal/informatics/workshops/ChIP-Seq/Slurm

The output files will appear on the filesystem at

    /n/regal/informatics/workshops/ChIP-Seq/Output
when the tutorials are run

## Prerequisites

* Ability to login to odyssey (ssh, openauth)
* Knowledge of basic command line commands (ls, cd, cat, less, nano, mkdir, |, >)
* Knowledge of basic script loops
* Knowledge of slurm batch scripts
* Ability to transfer files to/from the cluster


## Chip-Seq Summary
Chip-seq, or Chromatin Immunoprecipitation Sequencing, is used to probe DNA-protein interactions. First, DNA in whole cells is reversibly cross-linked to bound proteins with a cross-linking agent such as formaldehyde. Next, the DNA is isolated and sonicated to produce small fragments. Fragments containing the protein of interest are isolated via the addition of a protein-specific antibody and chromatin immunoprecipitation. Crosslinks are reversed and the isolated DNA is purified. Finally, strand-specific adapters are ligated and the library is amplified.

Strand-specific reads help to distinguish signal from noise, as true signal will consist of two peaks of equal magnitude, but on opposite strands. This is due to the fact that reads can originate from either end of DNA fragments, and DNA fragments are usually much longer than reads. This characteristic peak profile is shown in Figure 1.

<figure>
	<a class="img" href="/images/chipseq2.png">
    		<img class="img-responsive" src="/images/chipseq2.png"></img>
	</a>
    <figcaption>Figure 1</figcaption>
</figure>

## Outline of ChIP-Seq Sequencing Process

*  Design experiment and sequence
*  Pre-process the reads (trimmomatic)
*  Align the reads to a reference genome  (bowtie2)
*  Call peaks on the aligned reads using condition and control alignments (macs)
*  Investigate the content of the peaks for binding motifs (homer)
*  Search the peaks for a specific motif (homer)
*  Investigate the location of the peaks for enriched annotation and go terms (homer but probably not today)
	
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
These are essential to eliminate false positives as there is always some signal on open chromatin.


### Replicates
The good news here is that, unlike RNA-Seq, more than 2 replicates does not significantly increase the number of targets.


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
## Pre-processing and Alignment (Homework)
Due to time constraints we ask people to preprocess and align the raw sequence reads ahead of time.  This is meant as a refresher for basic file manipulation and slurm submission to reactivate your muscle memory for cluster use.    The walkthrough is available [here](http://informatics.fas.harvard.edu/chip-seq-workshop-alignment-homework/).

Pre-computed bam files are available in

    /n/regal/informatics/workshops/ChIP-Seq/Output/*.bam
If you want to use these files rather than your own you can link to them rather than copy them over.

1. If you haven't already done so log into the cluster and cd into the workshop directory
   <div class="codehilite"><pre>ssh login.rc.fas.harvard.edu
cd /n/regal/informatics/workshops/ChIP-Seq/Users</pre></div>

2. Make your own user directory if you haven't already and cd into it
   <div class="codehilite"><pre>mkdir &lt;username&gt;
cd &lt;username&gt;</pre></div>

3. Now you're going to do a fancy bash loop over the bam files and link to them so they look like they're in your directory
   <div class="codehilite"><pre>for i in /n/regal/informatics/workshops/ChIP-Seq/Output/*.bam* ; do
    ln -s $i .
done</pre></div>


Do a quick `ls -l` to check you've got the right links.  There should be bam files and bam.bai index files.

## MACs Peak Calling
We’re going to use MACS for peak calling which is a commonly used program and often comes out top in a comparison of methods.  The publication describing the algorithm is [here](http://www.genomebiology.com/content/9/9/R137).

Paramount in any peak calling method is the use of controls as reads are not distributed randomly across the genome and can generate many false positives.

In addition MACS empirically models the shift size of ChIP-Seq tags, and uses it to improve the spatial resolution of predicted binding sites. MACS also uses a dynamic Poisson distribution to effectively capture local biases in the genome.

## MACS Command Line
An example of a macs command line with commonly used options is below

    :::bash
    macs2 callpeak  \
    --name macs1
    --treatment file1.bam file2.bam \
    --control ctrl1.bam ctrl2.bam  \
    --outdir macs1 \
    --format BAM \
    --name macs1 \
    --pvalue 1e-5 \
    --gsize hs \
    --tsize 75 \
    --bdg \
    --wig

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

The only essential parameter is the `--treatment` one but I recommend you include the others (up to tsize) as a record of how you ran the program.

More parameters and their descriptions are at

[https://github.com/taoliu/MACS/blob/macs_v1/README.rst](https://github.com/taoliu/MACS/blob/macs_v1/README.rst)

## Slurm Submission
Although it's perfectly possible to run MACS interactively it's good practice to submit jobs through the job queueing system slurm.

Below is a template for a sample slurm script

    :::bash
    #!/bin/bash
    
    # Using $1 in the command means you can call this script with an input file
    #
    # sbatch thisfile.sh myinputfile.dat
    #
    # If you have multiple input files then you can use a bash for loop to submit them all
    #
    # for i in *.fastq ; do
    # sbatch thisfile.sh $i
    # done
    #
    #
    # Change the <nnnn> strings to something sensible before submitting
    #
    #SBATCH -J <jobname>
    #SBATCH -N 1                  # Ensure that all cores are on one machine
    #SBATCH -n <n>                # Use n cores for the job
    #SBATCH -t <n-nn:nn>          # Runtime in D-HH:MM
    #SBATCH -p <queuename>        # Partition to submit to
    #SBATCH --mem=<n>             # Memory pool for all cores in Mb (see also --mem-per-cpu)
    #SBATCH -o <outfile>.%A.out   # File to which STDOUT will be written (%A is replaced by the jobid)
    #SBATCH -e <outfile>.%A.err   # File to which STDERR will be written (%A is replaced by the jobid)
    #SBATCH --mail-type=<type>    # Type of email notification- BEGIN,END,FAIL,ALL
    #SBATCH --mail-user=<myemail@harvard.edu> # Email to which notifications will be sent

    # You can have multiple module loads here
    module load <mymodule>

    # The command can use the input parameter $1 here

    some_command_here $1 > $1.out
        
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

### Suggestions for MACS sbatch files
<table>
<tr><td>Running time</td><td>~1 hour per sample</td></tr>
<tr><td>Memory needs</td><td>~2G</td></tr>
<tr><td>Module</td><td><a href="https://portal.rc.fas.harvard.edu/apps/modules/macs2/2.1.0.20140616-fasrc01">macs2/2.1.0.20140616-fasrc01</a></td></tr>
<tr><td>Queue</td><td>serial_requeue usually schedules pretty quickly but jobs can get cancelled.</td></tr>
</table>

## Output files

The main output file is  NAME_peaks.narrowpeaks which is a  BED format file which contains the peak locations.   This will be the file you use as input for the motif analysis.

You can load it to UCSC genome browser or Affymetrix IGB software. The 5th column in this file is the -10*log10pvalue of peak region.

Other output formats are described at
[https://github.com/taoliu/MACS/blob/macs_v1/README.rst#output-files](https://github.com/taoliu/MACS/blob/macs_v1/README.rst#output-files)

## Exercise 1 - Peak Calling
1.  Construct an sbatch submission script for 1 treatment and control file and submit (this will run faster)
1.  Construct an sbatch submission script for all treatment and all control files and submit
1.  Use `squeue -u <username>` to check the jobs are submitted and either pending or running.   Look at the error files if if fails to run.
1.  Look at the output   (see /n/regal/informatics/workshops/ChIP-Seq/Output/macs2_* for pre-computed output)
    <ul><li>Download the peaks file, wig or bdg files and your bam files and indices to your laptop</li><li>Load up your bam files/wigs etc into IGV</li><li>Zoom and scroll through some peaks.   Look at some high scoring peaks and some low scoring peaks.  Is this what you expect?   Is there anything different (bar the obvious) between the high scoring and the low scoring peaks.</li><li>Extra:  The output files AF1.bed, AF2.bed from the HOXB7 paper are in the output directory.  Load these into IGV  and compare to your peaks.   How similar are they?   Can you write a command line to find the intersection between your results and the paper's results?</li><li>Extra:  Using the bedtools module  `bedtools2/2.25.0-fasrc01` use the bedtools intersect command to compare the AF1.bed and AF2.bed files to your results.</li></ul>
    

Hint:   To read the help info for a bedtools command do the following:

    :::bash
    bedtools intersect 2>&1 | more

Bedtools prints the help text to stderr (there are two output streams stdout and stderr).  By default `|more` will only redirect the stdout.   Using the `2>&1` will direct the stderr (2) output into the stdout stream (1)

## Motif Finding using Homer

[http://homer.salk.edu/homer/introduction/basics.html](http://homer.salk.edu/homer/introduction/basics.html)

Motif finding has many pitfalls.

*  The motifs are small and common in the genome,
*  The motifs are not always perfectly described in the database
*  There is a lot of overlap beween motifs (degeneracy)
*  Methods look for differential enrichment between a condition and control set of sequences.
*   Controls are hard to get right and need thought.

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
    
### Hands on</strong>

1. Construct and submit an sbatch file to find motifs on the MACS peaks for a motif of length 7 in a region of 50bp. Use multiple threads (16 is good),  1Gb memory and 30 mins run time.

2. Construct and submit an sbatch file to find motifs on the MACS peaks for motifs of lengths 7 in a region of 200bp.  Use multiple threads, 1Gb memory and 30 mins run time.<br/>
Either wait for the jobs to run or look at the pre-run output in

    `/n/regal/informatics/workshops/ChIP-Seq/Output/macs2_all/macs2_all_peaks.narrowPeak_homer_50bp_7mer` <br/><br/>
or for the 200bp region

    `/n/regal/informatics/workshops/ChIP-Seq/Output/macs2_all/macs2_all_peaks.narrowPeak_homer_100bp_7mer`
    
3. Transfer the output directory to your laptop and look at the html output.   Considering the protein used in this experiment are the results understandable?   Do the differences between the 50bp and 200bp region make sense.

4. How would you construct an sbatch file to search for motifs of multiple lengths?

## Using  Homer to search for a particular motif
The motif reported in the paper is GCCNGGC . How can we search for this motif in our data?

1. First we generate a motif file for our motif (no need for a slurm script for this)
    <div class="codehilite"><pre>seq2profile.pl &lt;consensus&gt; [# mismatches] [name] &gt; output.motif</pre></div>

2. Then we use that motif file to search our peaks file
<div class="codehilite"><pre>findMotifsGenome.pl  &lt;peakfile&gt;  \ 
    &lt;genome&gt; \ 
    &lt;outputdir&gt; \
    -size &lt;size&gt; \
    -p &lt;threads&gt; \
    -find &lt;motiffile&gt;
</pre></div>
  
### Hands on
We need to generate a motif file but we don't have a consensus.  Let's run with a fake consensus to generate a template.

    :::bash
    seq2profile.pl GCCAGGC  1 HOXB7 &gt; HOXB7.motif
    
Now we can edit the HOXB7.motif file to edit the probabilites to make the middle character 0.25,0.25,0.25,0.25

(Edit the A position to reflect the correct probabilities)

    >GCCAGGC    HOXB7    2.76827819473531
    0.001    0.001    0.997    0.001
    0.001    0.997    0.001    0.001
    0.001    0.997    0.001    0.001
    0.997    0.001    0.001    0.001
    0.001    0.001    0.997    0.001
    0.001    0.001    0.997    0.001
    0.001    0.997    0.001    0.001


Then we can run findMotifsGenome.pl with this file to find candidates

    :::bash
    findMotifsGenome.pl  macs2_peaks.narrowPeak hg19 macs2_homer_200_HOXB7 -size 200 -find HOXB7.motif  -p 16 > macs2_homer_200_HOXB7.out
    
### Hands on

1.  Run the `seq2profile.pl` command on the login node (no need for a slurm command for this).

2.  Edit the resulting file to have the correct probability distribution for the middle base.

3. Construct a slurm batch file for the findMotifsGenome.pl and submit to the cluster.  Use a run time of ~1hour and 1Gb RAM.

To save time the output is in :

    /n/regal/informatics/workshops/ChIP-Seq/Output/macs2_all/macs2_homer_200_HOXB7

How do your results compare to the results in the paper?     What do you think is different?