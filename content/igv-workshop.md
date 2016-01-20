Title: IGV Workshop
Date: 2015-01-01
Category: Tutorials
Author: Michele Clamp
Tags: IGV, Visualization
Summary: This tutorial provides useful tips for IGV analysis.

## Prerequisites and Installation.

Note:  The IGV user guide [http://www.broadinstitute.org/software/igv/UserGuide](http://www.broadinstitute.org/software/igv/UserGuide) is extremely good and has all this information and much more.  

*   Laptop (linux/osx/windows) with a bare minimum of 2Gb of RAM.   4Gb is better but there is no such thing as too much.
*   IGV installed.  You need to register your email address before you can get access by clicking on the ‘Register’ link and filling out the form.
*   A copy of the workshop data.  This will be distributed via usb drive at the start of the workshop.

<figure>
	<a class="img" href="/images/igv1.png">
    		<img class="img-responsive" src="/images/igv1.png"></img>
	</a>
    <figcaption></figcaption>
</figure>


<figure>
	<a class="img" href="/images/igv2.png">
    		<img class="img-responsive" src="/images/igv2.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Once you’ve filled out your details you’ll be taken to the download page.

<figure>
	<a class="img" href="/images/igv3.png">
    		<img class="img-responsive" src="/images/igv3.png"></img>
	</a>
    <figcaption></figcaption>
</figure>


### Other useful utilities - igvtools

Although not essential there is also a set of command line utilities igvtools that are useful.   If you download and unzip the file you can run them (assuming java is installed)

    :::bash
    unzip igvtools_2.3.60.zip
    ./IGVTools/igvtools

This will print out the help :

    :::bash
    [mclamp@Micheles-iMac:~/Downloads]$ ./IGVTools/igvtools
    
    Program: igvtools. IGV Version 2.3.60 (71)09/21/2015 06:08 PM
    
    Usage: igvtools [command] [options] [input file/dir] [other arguments]

    Command: version print the version number
         sort    sort an alignment file by start position.
         index   index an alignment file
         toTDF   convert an input file (cn, gct, wig) to tiled data format (tdf)
         count   compute coverage density for an alignment file
         formatexp  center, scale, and log2 normalize an expression file
         gui      Start the gui
         help <command></command>     display this help message, or help on a specific command
         See http://www.broadinstitute.org/software/igv/igvtools_commandline for more detailed help

    Error: No arguments provided
    Done


## The main functions

Before we fire up the program and load some data we’re going to go through the various parts of the main screen and what they do.IGV is a tool for viewing genomes, genome features and alignments to other sequences.    These features typically cover a range of bases and can either be

*   ‘plain’     i.e. just a chr:start-end with a name and a score (and maybe some other attributes).  Examples of these would be repeats,  CpG islands or SNPs.  File types for these include bed, gff, gff3, vcf, psl
*   ‘groups’     which are plain features but grouped into larger ones.  The obvious examples here are genes/transcripts which are groups of exons.
*   ‘alignments’ Typically for IGV these are short read alignments usually in bam format (although you can also use sam format).  The bam files need to be indexed before they can be read using a tool like samtools.

## The main screen

By default when you fire up IGV the most commonly used human genome version is loaded.  This is currently hg18.

<figure>
	<a class="img" href="/images/igv4.png">
    		<img class="img-responsive" src="/images/igv4.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

When the genome is first loaded all the chromosomes are shown in a line along with their names/numbers

<figure>
	<a class="img" href="/images/igv5.png">
    		<img class="img-responsive" src="/images/igv5.png"></img>
	</a>
    <figcaption></figcaption>
</figure>


Along the bottom is a density track for the associated annotation.  In the case of human these are refseq genes.

<figure>
	<a class="img" href="/images/igv6.png">
    		<img class="img-responsive" src="/images/igv6.png"></img>
	</a>
    <figcaption></figcaption>
</figure>


There are various buttons and controls along the top bar.   The most useful is on the far right which allows you to zoom in and out (although *not* in the whole genome view).    The other useful one is the home <img style="display: inline" src="/images/igv7.png"/> button which takes you back to this (whole genome) screen. Another part of the display to be aware of is the memory usage in the bottom right hand corner.  If you have large bam files and are displaying large numbers of reads this can quickly get used up and then IGV won't read any more reads.

## Narrowing the view down to a region

The whole genome view is not particularly useful and we usually want to narrow down to a region.  There are a number of ways to do this.

a) To focus on a single chromosome you can either click on one of the chromosome identifiers or use the dropdown box next to the genome dropdown.

<figure>
	<a class="img" href="/images/igv8.png">
    		<img class="img-responsive" src="/images/igv8.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

The view changes in a number of ways.

<ul>
<li>The whole genome view changes to a karyotype view of the chromosome along with size information.

<figure>
	<a class="img" href="/images/igv9.png">
    		<img class="img-responsive" src="/images/igv9.png"></img>
	</a>
    <figcaption></figcaption>
</figure>
</li>
<li>The gene track now shows genes rather than a coverage plot.

<figure>
	<a class="img" href="/images/igv10.png">
    		<img class="img-responsive" src="/images/igv10.png"></img>
	</a>
    <figcaption></figcaption>
</figure>
</li>
<li>The zoom controls are now active and you can zoom in and out using the plus and minus buttons.

<figure>
	<a class="img" href="/images/igv11.png">
    		<img class="img-responsive" src="/images/igv11.png"></img>
	</a>
    <figcaption></figcaption>
</figure>
</li>
</ul>

b) To focus on a known region of a chromosome

You can type a region string into the text field next to the chromosome dropdown.  For instance you can type `chr10:1100000-10000000`

<figure>
	<a class="img" href="/images/igv12.png">
    		<img class="img-responsive" src="/images/igv12.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Note that a red box appears on the karyotype showing you where you are on the chromosome and the blue box on the zoom controls shows you how zoomed in you are. You don’t have to put in a range - you can just type a base chr10:1100000   

<figure>
	<a class="img" href="/images/igv13.png">
    		<img class="img-responsive" src="/images/igv13.png"></img>
	</a>
    <figcaption></figcaption>
</figure>


Note this only shows 41bp which is a little small for most purposes.Finally you can type in a gene name or identifier that is in the annotation track e.g. ALX1.  As you type you’ll get a list of possible options. 

<figure>
	<a class="img" href="/images/igv14.png">
    		<img class="img-responsive" src="/images/igv14.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

And when you select the gene the region display changes to the base range. 

<figure>
	<a class="img" href="/images/igv15.png">
    		<img class="img-responsive" src="/images/igv15.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

### Other navigation methods

*   Clicking on the karyotype will center the display at that point.
*   Double clicking a feature will center the display on that feature.
*   The plus and minus keys zoom in and out.
*   Shift click and option click zoom in and out.

### Scrolling

There are numerous ways to scroll through the display

*   You can drag the panels left and right using the mouse
*   You can move a short way using the left and right arrows.
*   You can move a whole screen left and right using the home and end keys.

## Loading IGV Provided Genomes

To load up a provided genome (or to search to see if it exists) go to the Genome dropdown in the top left and click 'More...' 

<figure>
	<a class="img" href="/images/igv16.png">
    		<img class="img-responsive" src="/images/igv16.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

This will bring up a window where you can scroll through a list of provided genomes or use the search field to search for what you need. 

<figure>
	<a class="img" href="/images/igv17.png">
    		<img class="img-responsive" src="/images/igv17.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

If you find it IGV will download both the karyotype and the annotation files and display them for you.  

<div class="exercises">
Exercises
* Load up the Zebrafish Zv9 genome.
* Which chromosome is the ALX4b gene on?
* Which genes are upstream and downstream of the ALX4b gene
</div>

## User-provided Genomes

If your genome isn’t in the list you can still use IGV if you have a reference genome fasta file with one sequence per chromosome.   You can also provide an annotation set.Finding your genome and annotation - coming soon!

<figure>
	<a class="img" href="/images/igv18.png">
    		<img class="img-responsive" src="/images/igv18.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

This will just add the fasta filename to your genomes list.

<figure>
	<a class="img" href="/images/igv19.png">
    		<img class="img-responsive" src="/images/igv19.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Alternatively if you have an annotation file as well you can easily create a custom genome using the Genomes->Create .genome file option 

<figure>
	<a class="img" href="/images/igv20.png">
    		<img class="img-responsive" src="/images/igv20.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Fill out the fields and select the genome fasta and the annotation gff or gtf file and you get a named genome and a gene track as we did for the human genome.

<figure>
	<a class="img" href="/images/igv21.png">
    		<img class="img-responsive" src="/images/igv21.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

<figure>
	<a class="img" href="/images/igv22.png">
    		<img class="img-responsive" src="/images/igv22.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

<div class="exercises">
Exercises
<ul>
<li>Create a .genome file using the taeGut3.fa and taeGut3.ensembl82.gff files in your IGV_Workshop/Data/taeGut3 directory.</li>
<li>Which chromosome is the ALX1 gene on?</li>
<li>Which genes are upstream and downstream of the ALX1 gene?</li>
</ul>
</div>

## Loading tracks/features and alignments

So now we’re going to look at some data from a ChIP-Seq experiment.  For this we have the following

*  genome (hg19)
*  aligned read data (.bam files) from a sample and a control.
*  Called peak regions.

First we select hg19 from the genome menu.  If it isn’t there select the ‘More…’ option and search for it. 

<figure>
	<a class="img" href="/images/igv23.png">
    		<img class="img-responsive" src="/images/igv23.png"></img>
	</a>
    <figcaption></figcaption>
</figure>


Note we have chromosomes and a gene track automatically displayed.  Now we’re going to load up 2 bam files for the sample and the control from your IGV_Workshop/Data/ChIP-Seq directory.

* sample file is `SRR866428_1.fastq.trimmed.fastq.bam` 
* control file is `SRR866432_1.fastq.trimmed.fastq.bam`

Note that the index files SRR866428_1.fastq.trimmed.fastq.bam.bai etc need to be in the same directory for IGV to be able to read the bam files. If they’re not there you need to use samtools to create them on the cluster using :

    :::bash
    module load samtools/0.1.19-fasrc01
    samtools index myfile.bam

<figure>
	<a class="img" href="/images/igv24.png">
    		<img class="img-responsive" src="/images/igv24.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

If this is successful you’ll see two new tracks but no alignments as at this zoomlevel there are too many to show.  We have to zoom in to see the reads.

<div class="exercises">
Exercise
<ul>
<li>Load the hg19 genome into IGV</li>
<li>Load the two bam files</li>
<li>At what zoom level (how many kb) do the reads appear?</li>
</ul>
</div>

## Read information and color schemes

Let’s look at the region chr1:91,334,961-91,335,725

<figure>
	<a class="img" href="/images/igv25.png">
    		<img class="img-responsive" src="/images/igv25.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

*   We can see that most of the reads are colored gray.  This means that their mapping quality is 1.   If a read is badly mapped it will be transparent.
*   The colored vertical lines show base mismatches to the reference.
*   Detailed information on the read and each base is shown if you mouse over a read.

<figure>
	<a class="img" href="/images/igv26.png">
    		<img class="img-responsive" src="/images/igv26.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Many more options for displaying reads can be seen if you right click in the track window. 

<figure>
	<a class="img" href="/images/igv27.png">
    		<img class="img-responsive" src="/images/igv27.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Different color schemes can be selected from the 'Color alignments by’ sub-menu.  For single ended chip-seq data coloring by strand is a good option. 

<figure>
	<a class="img" href="/images/igv28.png">
    		<img class="img-responsive" src="/images/igv28.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Other useful options for single end reads are switching between Collapsed/Expanded/Squished to change how the reads are displayed.  We’ll investigate more of the options for paired-end data in the next section.

## Loading feature files

For ChIP-Seq data we generally have a set of called peak regions.   We have a set of these in bed format in the AF1.bed file.  We can load these up and display them from the File->Load from File… menu and then select the AF1.bed file from you IGV_Workshop/Data/ChIP-Seq directory 

<figure>
	<a class="img" href="/images/igv29.png">
    		<img class="img-responsive" src="/images/igv29.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

<figure>
	<a class="img" href="/images/igv30.png">
    		<img class="img-responsive" src="/images/igv30.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

This should display a new track at the bottom of the panel. 

<figure>
	<a class="img" href="/images/igv31.png">
    		<img class="img-responsive" src="/images/igv31.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

There aren’t that many peaks in this file and we’d like to show them all at whatever zoom level we’re at.   We can change this by right clicking on the AF1 track and selecting ‘Set Feature Visibility Window….' 

<figure>
	<a class="img" href="/images/igv32.png">
    		<img class="img-responsive" src="/images/igv32.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Enter 0 to always display all data. 

<figure>
	<a class="img" href="/images/igv33.png">
    		<img class="img-responsive" src="/images/igv33.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

### Moving from feature to feature

As these peaks aren’t very dense scrolling is somewhat tedious to view them.    You can move from feature to feature by selecting the track name (AF1) and using ctrl-F and ctrl-B to move forward to the next feature or back to the previous feature. 

<figure>
	<a class="img" href="/images/igv34.png">
    		<img class="img-responsive" src="/images/igv34.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

<div class="exercises">
Exercise
<ul>
<li>Load up the AF1.bed file</li>
<li>Set the feature visibility to 0</li>
<li>Look at the read patterns for a few peaks.  Change the color scheme for the reads to by strand.   Do the peaks look as though they have 50/50 forward/reverse reads.</li>
</ul>
</div>

### Having a list of regions

We can also load up the peaks as a scrollable list as an alternative way of navigating through them.   To do this go to  Regions->Import Regions and select the AF1.bed file. 

<figure>
	<a class="img" href="/images/igv35.png">
    		<img class="img-responsive" src="/images/igv35.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

After that go back to the Region menu and choose Region Navigator…   You should now see a window with all your peaks, their positions and names.

<figure>
	<a class="img" href="/images/igv36.png">
    		<img class="img-responsive" src="/images/igv36.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

  From here you can either select one at a time to view or you can select multiples and they will display side by side.

<figure>
	<a class="img" href="/images/igv37.png">
    		<img class="img-responsive" src="/images/igv37.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

<figure>
	<a class="img" href="/images/igv38.png">
    		<img class="img-responsive" src="/images/igv38.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

To go back to a single region view double click on the id string (chr19 in this case) at the top of the window.

<div class="exercises">
Exercises
<ul>
<li>Load up the AF1.bed file as a region list.</li>
<li>Select some peaks and display in a multi region panel.</li>
<li>Go back to a single pane view.</li>
</ul>
</div>

## Paired-end RNA-Seq data

For the next section we’re going to look at paired end data from an RNA-Seq run.   The data we have is

*   A Zebrafinch reference genome with annotation (you created this in the first part of the workshop)
*   A de novo (Trinity) assembly of RNA-Seq data from a close but non-identical genome (GenomeX) aligned to the Zebrafinch genome
*   A tophat alignment of the GenomeX reads to the Zebrafinch genome.

Our aim is to make a qualitative decision as to whether the de novo Trinity assembly is better than the tophat2 alignment for closely related genomes (~95% identical).

First we switch to the Zebrafinch assembly we created previously using the dropdown in the top left of our IGV window. 

<figure>
	<a class="img" href="/images/igv39.png">
    		<img class="img-responsive" src="/images/igv39.png"></img>
	</a>
    <figcaption></figcaption>
</figure>


### Switch on splice junction display.

IGV has a nice feature where it can show reads spanning introns more clearly than just looking at the read alignments.   This isn’t on by default so we need to go into the preferences panel to turn it on.   Choose the View->Preferences menu option and the preference panel will be displayed.  From here Choose the ‘Alignment’ tab. 

<figure>
	<a class="img" href="/images/igv40.png">
    		<img class="img-responsive" src="/images/igv40.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

<figure>
	<a class="img" href="/images/igv41.png">
    		<img class="img-responsive" src="/images/igv41.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

In the second section from the bottom select the ‘show junction track’ and ‘show flanking regions’ and click ok.  We now want to load our tophat alignment and our de novo assembly transcript alignments using the File->Load from File… menu

*   Tophat file is a bam file TaeGut3/tophat2_taeGut3/accepted_hits.bam
*   Trinity alignment file is a psl file from blat TaeGut3/Trinity.blat.tophit.sort.psl

<div class="exercises">
Exercises
<ul>
<li>Load the ZebraFinch genome you created earlier.</li>
<li>Set the preferences to show junctions and flanking regions</li>
<li>Load the TaeGut3/tophat2_taeGut3/accepted_hits.bam file</li>
<li>Load the TaeGut3/Trinity.blat.tophit.sort.psl file</li> </li>
<li>Go to the ALX1 gene</li>
<li>Change the various track displays (expand/collapse, track height etc) so individual transcripts can be seen.</li>
</ul>
</div>

You should end up with something like

<figure>
	<a class="img" href="/images/igv42.png">
    		<img class="img-responsive" src="/images/igv42.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

### Changing the display of paired end data

As with the single end reads the default display is by mapping quality (transparent is 0) and highlighting mismatched bases.  As this is paired end data if a read straddles an exon a line is drawing joining the pieces.

Intron spanning regions are also shown in the Junctions track.   the height of the spanning curve indicates how many reads span that intron and the color denotes strand (red = forward, blue = reverse).

Again if you mouse over a read you get a popup with read and alignment info.   This time you also get information on the mate pair.

<figure>
	<a class="img" href="/images/igv43.png">
    		<img class="img-responsive" src="/images/igv43.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

The junction curves also show info if you mouse over them too. 

<figure>
	<a class="img" href="/images/igv44.png">
    		<img class="img-responsive" src="/images/igv44.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

If you want to find the mate pair of a read you can option-click (ctrl-click on windows) and it will be highlighted in the same color.  

If you now right click the menu brought up now has extra options that were greyed out previously. 

<figure>
	<a class="img" href="/images/igv45.png">
    		<img class="img-responsive" src="/images/igv45.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

The first one we’ll try is ‘view as pairs’.  This reorganizes the display so mate pairs are shown one after the other.  The mouseover also changes to display information about both reads in the pair. 

<figure>
	<a class="img" href="/images/igv46.png">
    		<img class="img-responsive" src="/images/igv46.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Changing the color scheme can give some insight into alignment quality.  Choosing the color by ‘insert size and pair orientation’ will color reads as

*   Blue is for inserts that are smaller than expected
*   Red is for inserts that are larger than expected.
*   Inter-chromosomal rearrangements are color-coded by chromosome.
*   Shades of green, teal, and dark blue show structural events of inversions, duplications, and translocations.


If the mate pair is somewhere else in the genome pick ‘view mate in split screen’ from the right click popup menu. (Ran) 

<figure>
	<a class="img" href="/images/igv47.png">
    		<img class="img-responsive" src="/images/igv47.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

<figure>
	<a class="img" href="/images/igv48.png">
    		<img class="img-responsive" src="/images/igv48.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

<div class="exercises">
Exercises
<ul>
<li>Switch your display to ‘view as pairs’ and color by insert size and orientation.</li>
<li>Scroll through a chromosoem gene by gene (select the Gene track name and ctrl-f ctrl-b to go back and forth) to answer the following:</li>
<li> Can you find a convincing structural rearrangment (many reads mapping to a different chromosome/wrong orientation)</li>
<li> Are the blat alignments a better annotation (compared to the Zebrafinch) than the splice junctions from tophat?</li>
</ul>
</div>
