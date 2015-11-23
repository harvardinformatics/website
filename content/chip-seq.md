Title: Chip-Seq
Date: 2015-01-02 11:00
Category: Tutorials
Tags: Chip-Seq
Summary: Chip-seq, or Chromatin Immunoprecipitation Sequencing, is used to probe DNA-protein interactions. This tutorial illustrates an example Chip-seq analysis. After reads are preprocessed and aligned, peaks are called and motifs identified. 

Chip-seq, or Chromatin Immunoprecipitation Sequencing, is used to probe DNA-protein interactions. First, DNA in whole cells is reversibly cross-linked to bound proteins with a cross-linking agent such as formaldehyde. Next, the DNA is isolated and sonicated to produce small fragments. Fragments containing the protein of interest are isolated via the addition of a protein-specific antibody and chromatin immunoprecipitation. Crosslinks are reversed and the isolated DNA is purified. Finally, strand-specific adapters are ligated and the library is amplified.

Strand-specific reads help to distinguish signal from noise, as true signal will consist of two peaks of equal magnitude, but on opposite strands. This is due to the fact that reads can originate from either end of DNA fragments, and DNA fragments are usually much longer than reads. This characteristic peak profile is shown in Figure 1.

[![chipseq2](http://informatics.fas.harvard.edu/wp-content/uploads/2014/06/chipseq2.png)](http://informatics.fas.harvard.edu/wp-content/uploads/2014/06/chipseq2.png)

Figure 1: Characteristic profile of chip-seq reads. (Park 2009 Nat. Reviews [PMID: 19736561](http://sfx.hul.harvard.edu/hvd?url_ver=Z39.88-2004&rft_val_fmt=info:ofi/fmt:kev:mtx:journal&__char_set=utf8&rft_id=info:pmid/19736561&rfr_id=info:sid/libx%3Ahul.harvard&rft.genre=article))

As the cost of sequencing has decreased, it has become more feasible to include replicate and control samples in Chip-seq experiments. Typically, control samples are prepared in the same manner as experimental samples, except that no immunoprecipitation is performed. Alternatively, immunoprecipitation may be performed without antibody, or with an antibody not known to interact with DNA, such as Immunoglobulin G. By subjecting the control sample to the same processing steps, but without selecting for a specific protein, control samples provide a read count profiles that represent local chromatin structure, genome copy variation and sequence bias. They also help to filter out peaks in read counts due to read mis-alignment.

The following scripts illustrate an example Chip-seq analysis. After reads are preprocessed and aligned, peaks are called and motifs identified. The scripts below show typical parameter settings as well as the CPU, memory and time requirements of each step.

### Preprocessing

The first step in Chip-seq analysis is to filter low-quality reads and trim adapter sequences from read ends. The following script uses Trimmomatic to remove Illumina single-end TruSeq adapters and trim low-quality bases. In addition, it enforces a minimum quality score over a sliding window and requires a minimum read length:

<pre>#!/bin/bash
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=16                    # Number of cores
#SBATCH --time=2:00:00                 # Runtime HH:MM:SS
#SBATCH --partition=serial_requeue     # Partition to submit to
#SBATCH --mem-per-cpu=2000             # Memory in MB     

module load centos6/Trimmomatic-0.30
java -jar $TRIMMOMATIC/trimmomatic-0.30.jar SE \
     -threads 16 \
     -phred33 \
     Sample1.fastq Sapmple1_pp.fastq \
     ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq3-SE.fa:2:30:10 \
     LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36</pre>

To run this job on the cluster, save the script in a file called preprocess.sh, and run the command:

<pre>sbatch preprocess.sh</pre>

For more information on using SLURM commands to run jobs, see this [documentation](https://rc.fas.harvard.edu/resources/running-jobs/).

### Mapping to a Reference Genome

After the reads are processed, they can be mapped to a reference genome with BWA:

<pre>#!/bin/bash
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=16                    # Number of cores
#SBATCH --time=2:00:00                 # Runtime HH:MM:SS
#SBATCH --partition=serial_requeue     # Partition to submit to
#SBATCH --mem=5000                     # Memory in MB

module load centos6/bwa-0.7.7
bwa aln -t 16 hg19.fa Sample1.fastq > Sample1.sai
bwa samse hg19.fa Sample1.sai Sample1.fastq > Sample1.sam</pre>

###  Peak Calling

The simplest way to call peaks would be to subtract control from experimental read counts at each location (if a control sample was present) and then compute a sliding window sum of read counts in the experimental sample. A more robust technique would examine read strandedness for the characteristic two-peak pattern mentioned above, using an estimation of the expected distance between peaks. In addition, the control sample(s) would be used to derive a local background model representing local chromatin structure, genome copy variation and other sources of background signal such as sequence bias, read mis-alignment and non-specific pull-down. We can perform these steps using the MACS algorithm, as shown below:

<pre>#!/bin/bash
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of cores
#SBATCH --time=2:00:00                 # Runtime HH:MM:SS
#SBATCH --partition=serial_requeue     # Partition to submit to
#SBATCH --mem=4000                     # Memory in MB

module load centos6/MACS-2.0.10_python-2.7.3
macs2 callpeak --treatment Sample1.sam Sample2.sam \
               --control Control1.sam Control2.sam \
               --name HypoxiaVsInput \
               --format SAM \
               --gsize hs \
               --tsize 50</pre>

The above example, there are two experimental replicates, and two replicate control samples. The read length, or tag size, is specified with the tsize parameter. If using a common reference genome, the size of the genome can be specified with a shortcut such as “hs,” as shown above. Otherwise, the number of bases can be directly specified.

### Peak Visualization

MACS generates a list of peaks in [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) format, which can be visualized in [IGV](https://www.broadinstitute.org/igv/home) along with the aligned reads from which the peaks were generated:

<div id="attachment_939" style="width: 784px" class="wp-caption alignnone">[![chipseq3](http://informatics.fas.harvard.edu/wp-content/uploads/2014/06/chipseq3.png)](http://informatics.fas.harvard.edu/wp-content/uploads/2014/06/chipseq3.png)

Figure 2: Aligned reads and called peaks in IGV. (Feng 2012 Nat. Prot [PMID: 22936215](http://sfx.hul.harvard.edu/hvd?url_ver=Z39.88-2004&rft_val_fmt=info:ofi/fmt:kev:mtx:journal&__char_set=utf8&rft_id=info:pmid/22936215&rfr_id=info:sid/libx%3Ahul.harvard&rft.genre=article)

</div>

### Motif Discovery

Once peaks have been identified, we can look for motifs that are enriched in peaks compared to the background genome sequence. he [MEME](http://meme.nbcr.net/meme/) suite allows for novel motif discovery as well as comparison of found motifs to a database of known motifs. Because of MEME’s database integration, and typically small input file size of peak sequences, motif discovery is often run through their web [submission form](http://meme.nbcr.net/meme/cgi-bin/meme.cgi).

[![chipseq4](http://informatics.fas.harvard.edu/wp-content/uploads/2014/06/chipseq4.png)](http://informatics.fas.harvard.edu/wp-content/uploads/2014/06/chipseq4.png)
