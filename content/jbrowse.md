Title: JBrowse on the FASRC Cluster
Date: 2019-06-19
Modified: 2021-06-23
Author: Nathan Weeks
Category: Software
Tags: Genome Annotation, JBrowse, Open OnDemand
Summary: How to launch JBrowse on the FASRC Cluster with Open OnDemand

[TOC]

## Introduction

[JBrowse](https://jbrowse.org/) is a web-based genome browser for visualizing genomic features in common file formats, such as variants (VCF), genes (GFF3, BigBed) and gene expression (BigWig), and sequence alignments (BAM, CRAM, and GFF3).

Both [JBrowse 1](https://jbrowse.org/jbrowse1.html) and [JBrowse 2](https://jbrowse.org/jb2/) can be launched from the FAS RC VDI (see [Running JBrowse on the FASRC Cluster using Open OnDemand](#running-jbrowse-on-the-fasrc-cluster-using-open-ondemand)).
JBrowse 2 is generally recommended for new use.

## JBrowse 2 configuration file (config.json)

The [jbrowse CLI](https://jbrowse.org/jb2/docs/cli/) can be used to generate the JBrowse 2 [config.json](https://jbrowse.org/jb2/docs/config_guide#intro-to-the-configjson).
A [Singularity](https://docs.rc.fas.harvard.edu/kb/singularity-on-the-cluster/) image containing the jbrowse/cli npm module is available on the FAS RC cluster, and used in the example below.

The following sequence of commands initializes a config.json with a set of reference sequences, adds two tracks of files that refer to the reference sequence IDs, and sets the default session to a linear view of the chromosomes (note all files are assumed to be formatted appropriately; see [JBrowse data file formats](#jbrowse-data-file-formats)):

    :::bash
    jbrowse_cli=/n/singularity_images/informatics/jbrowse/jbrowse2_v1.3.0.sif
    singularity exec ${jbrowse_cli} jbrowse add-assembly reference.fa.gz --load inPlace
    singularity exec ${jbrowse_cli} jbrowse add-track genes.sorted.gff3.gz --load inPlace
    singularity exec ${jbrowse_cli} jbrowse add-track variants.sorted.vcf.gz --load inPlace
    singularity exec ${jbrowse_cli} jbrowse set-default-session -v LinearGenomeView

Note that all data files must exist in the same directory, or subdirectory of, the directory with the config.json file (note that hard links to files located elsewhere on the same file system are possible).

## JBrowse 1 configuration file (tracks.conf)

The following is an example of a tracks.conf [minimal configuration](https://jbrowse.org/docs/minimal.html) that lets JBrowse 1 infer track types and defaults from file name suffixes (see entries under "Configuring Tracks" section at the aforementioned minimal configuration link for additional configuration options).

```
[GENERAL]
refSeqs=reference.fa.gz.fai
# reference sequence compressed with bgzip and indexed with "samtools faidx"

[tracks.reference]
urlTemplate=reference.fa.gz
# file specified as the refSeqs value without the ".fai" suffix

[tracks.alignments]
urlTemplate=alignments.sorted.bam
# my.sorted.bam.bai (created by "samtools index") is assumed to exist
## baiUrlTemplate=alignments.sorted.bai
# custom name if corresponding .bai file isn't suffixed with .bam.bai

[tracks.variants]
urlTemplate=variants.sorted.vcf.gz
# variants.sorted.vcf.gz.tbi assumed to exist

[tracks.genes]
urlTemplate=genes.sorted.gff3.gz
# genes.sorted.gff3.gz.tbi assumed to exist

[tracks.expression]
urlTemplate=expression.bw

[tracks.peaks]
urlTemplate=peaks.bed.gz
# peaks.bed.gz.tbi assumed to exist
```

`tracks.conf` and all data files listed therein are assumed to be in the same directory (typically named "data" and referred to as the JBrowse data directory; note that data files may be in subdirectories); e.g.:

```
$ ls data/
alignments.sorted.bam
alignments.sorted.bam.bai
config.json     # JBrowse 2
expression.bw
genes.sorted.gff3.gz
genes.sorted.gff3.gz.tbi
reference.fa.gz
reference.fa.gz.fai
reference.fa.gz.gzi
tracks.conf     # JBrowse 1
variants.sorted.vcf.gz
variants.sorted.vcf.gz.tbi
```

## JBrowse data file formats

For optimal performance on FASRC [networked, high-performance storage](https://www.rc.fas.harvard.edu/resources/odyssey-storage/#Networked_High-performance_Shared_Scratch_Storage),  data backing JBrowse tracks should be stored in compressed, indexed file formats that represent data for a given track in a single file (plus one or two smaller associated index files).
Examples are below.

In particular, the JBrowse 1 [flatfile-to-json.pl](https://jbrowse.org/docs/flatfile-to-json.pl.html) script should *not* be used, as it converts the input GFF3/GenBank/BED file into numerous smaller files that degrade performance on network/parallel file systems.

### Software Environment Setup

For the commands listed below, `bgzip` and `tabix` are provided on the FASRC cluster by the `htslib` environment module, while `samtools` and `bcftools` are provided by their own environment modules:
    :::bash
    module load bcftools htslib samtools

### GFF3

A [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) file should first be sorted by seqid (column 1) and start (column 4), then [bgzip](https://www.htslib.org/doc/bgzip.html)-compressed, and finally [tabix](https://www.htslib.org/doc/tabix.html)-indexed.

```
awk -F '\t' 'NF == 9' genes.gff3 | sort -k 1,1 -k 4n,4n | bgzip > genes.sorted.gff3.gz
tabix -p gff genes.sorted.gff3.gz
```

The above command generates a GFF3 lacking directives/comments found in the original, but they not necessary for viewing in JBrowse.

For very large GFF3 files, (GNU coreutils) sort options such as `--parallel=8 --buffer-size=1G --stable --compress-program=gzip` can be used to improve sorting performance, while bgzip `--compress-level=9` option may be used for extra compression, and `--threads=...` option can be used to specify additional threads for faster compression.

### BED

A number of [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file formats can be displayed after sorting by chrom (column 1) then chromStart (column 2), compressing with [bgzip](https://www.htslib.org/doc/bgzip.html), and indexing with [tabix](https://www.htslib.org/doc/tabix.html).

```
sort -k 1,1 -k 2n,2n peaks.bed | bgzip > peaks.sorted.bed.gz
tabix -p bed peaks.sorted.bed.gz
```

### BAM / CRAM

BAM and CRAM files must be sorted with `samtools sort` and indexed with `samtools index`.
samtools is available in the [Bioconda samtools package](https://bioconda.github.io/recipes/samtools/README.html).

```
samtools sort -o alignments.sorted.bam my.bam
samtools index alignments.sorted.bam
```

For large BAM files, use SAMtools version >= 1.6; also, see the `--threads` option to set the number of threads for sorting, the `-m` option to specify the amount of memory per thread, and the `-l` option to set compression level for the resulting BAM file.

### VCF

VCF files should be sorted with [bcftools](https://www.htslib.org/doc/bcftools.html#sort) (available in the [Bioconda bcftools package](https://bioconda.github.io/recipes/bcftools/README.html)), and indexed with [tabix](https://www.htslib.org/doc/tabix.html) (available in the [Bioconda htslib package](https://bioconda.github.io/recipes/htslib/README.html)).

```
bcftools sort --output-type=z --output-file=variants.sorted.vcf.gz variants.vcf
tabix -p vcf variants.sorted.vcf.gz
```

For large VCF files (above ~768M), the bcftools `--max-mem` option may be used to allocate extra memory for sorting VCF records, avoiding the use of the file system for an [external merge sort](https://en.wikipedia.org/wiki/External_sorting#External_merge_sort).

## Making feature names searchable in JBrowse

The JBrowse [generate-names.pl](https://jbrowse.org/docs/generate-names.pl.html) script (available in the [Bioconda jbrowse package](https://bioconda.github.io/recipes/jbrowse/README.html)) generates an index of feature names for searching and autocompletion when typed into the JBrowse search box, or the *View > search for features* menu option. This script may be invoked as demonstrated below (note the $JBROWSE_SOURCE_DIR environment variable is set on activation of the conda environment that jbrowse has been installed in, and the `--hashBits 4` option is an optimization to reduce the number of index files generated).

```
# Execute from the JBrowse data/ directory.
# "genes" and "variants" correspond to track identifiers in tracks.conf; i.e. [tracks.genes] and [tracks.variants]
${JBROWSE_SOURCE_DIR}/bin/generate-names.pl --tracks genes,variants --hashBits 4 --compress --out .
```

## Running JBrowse on the FASRC Cluster using Open OnDemand

A JBrowse instance can be launched on the FASRC cluster using the [Open OnDemand web portal](https://docs.rc.fas.harvard.edu/kb/virtual-desktop/) ([https://vdi.rc.fas.harvard.edu/](https://vdi.rc.fas.harvard.edu/)).
From the menu bar, select *Interactive Apps > JBrowse*.
The web form will be used to submit a SLURM job that runs on a compute node (see [https://docs.rc.fas.harvard.edu/kb/running-jobs/#Slurm_partitions](https://docs.rc.fas.harvard.edu/kb/running-jobs/#Slurm_partitions) for a list of some the commonly-available partitions; lab-specific partitions may also be specified, as well as a comma-separated list of partitions that SLURM may consider when scheduling the job).
The custom JBrowse Singularity container image runs the [lighttpd](https://www.lighttpd.net/) web server, which requires only a single processor core and minimal memory (256 MB is requested for the purpose SLURM job scheduling).

In the "path of a JBrowse data directory" textbox, enter the absolute path to the desired JBrowse data/ directory containing the config.json (if using JBrowse 2) or tracks.conf (if using JBrowse 1) created previously, then click "Launch".
Once the JBrowse job has been scheduled and the JBrowse instance started on a compute node, a "Connect to JBrowse" button will appear; click to open JBrowse in a new tab.

## Support

For assistance with launching JBrowse on the FASRC cluster, please submit a [help ticket](https://portal.rc.fas.harvard.edu/rcrt/submit_ticket).

JBrowse bug reports/feature requests may be submitted to the [JBrowse 1 issue tracker](https://github.com/GMOD/jbrowse/issues) or [JBrowse 2 issue tracker](https://github.com/GMOD/jbrowse-components) on GitHub.

For general JBrowse questions, see the [JBrowse Gitter chatroom](https://gitter.im/GMOD/jbrowse).
