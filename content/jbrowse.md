Title: JBrowse on Odyssey
Date: 2019-06-12
Author: Nathan Weeks
Category: Tutorials
Tags: Odyssey, Genome Annotation
Summary: How to use Open OnDemand to launch JBrowse on Odyssey

## Introduction

[JBrowse](https://jbrowse.org/) is a web-based genome browser for visualizing genomic features in common file formats, such as variants (VCF), genes (GFF3, BigBed) and gene expression (BigWig), and sequence alignments (BAM, CRAM, and GFF3).

## JBrowse configuration file (tracks.conf)

The following is an example of a tracks.conf [minimal configuration](https://jbrowse.org/docs/minimal.html) that lets JBrowse infer track types and defaults from file name suffixes (see entries under "Configuring Tracks" section at the aforementioned minimal configuration link for additional configuration options).

```
[GENERAL]
refSeqs=reference.fa.gz.fai            # reference sequence compressed with bgzip and indexed with "samtools faidx"
[tracks.reference]
urlTemplate=reference.fa.gz            # file specified as the refSeqs value without the ".fai" suffix
[tracks.alignments]
urlTemplate=alignments.sorted.bam      # my.sorted.bam.bai (created by tabix) is assumed to exist
# baiUrlTemplate=alignments.sorted.bai # custom name if corresponding .bai file isn't suffixed with .bam.bai
[tracks.variants]
urlTemplate=variants.sorted.vcf.gz     # variants.sorted.vcf.gz.tbi assumed to exist
[tracks.genes]
urlTemplate=genes.sorted.gff3.gz       # genes.sorted.gff3.gz.tbi assumed to exist
[tracks.expression]
urlTemplate=expression.bw
```

All files listed as urlTemplate values are assumed to be in the same directory (a directory typically named "data" and referred to as the JBrowse data directory) as tracks.conf (note that symbolic links to files accessible from Odyssey compute nodes may be used instead); e.g.:

```
$ ls data/
alignments.sorted.bam
alignments.sorted.bam.bai
expression.bw
genes.sorted.gff3.gz
genes.sorted.gff3.gz.tbi
reference.fa.gz
reference.fa.gz.fai
reference.fa.gz.gzi
tracks.conf
variants.sorted.vcf.gz
variants.sorted.vcf.gz.tbi
```

## JBrowse data file formats

For optimal performance on Odyssey's [networked and parallel file systems](https://www.rc.fas.harvard.edu/resources/odyssey-storage/), data backing JBrowse tracks should be stored in compressed, indexed file formats that represent data for a given track in a single file (plus one or two smaller associated index files).
Examples are below.

In particular, the JBrowse [flatfile-to-json.pl](https://jbrowse.org/docs/flatfile-to-json.pl.html) script should *not* be used on Odyssey, as it converts the input GFF3/GenBank/BED file into numerous smaller files that degrade performance on network/parallel file systems.

### GFF3

A [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) file should first be sorted by seqid (column 1) and start (column 4), then [bgzip](https://www.htslib.org/doc/bgzip.html)-compressed, and finally [tabix](https://www.htslib.org/doc/tabix.html)-indexed.

```
awk -F '\t' 'NF == 9' genes.gff3 | sort -k 1,1 -k 4n,4n | bgzip > genes.sorted.gff3.gz
tabix -p gff genes.sorted.gff3.gz
```

The above command generates a GFF3 lacking directives/comments found in the original, but they are not necessary for viewing in JBrowse.

For very large GFF3 files, (GNU coreutils) sort options such as `--parallel=8 --buffer-size=1G --stable --compress-program=gzip` can be used to improve sorting performance, while bgzip `--compress-level=9` option may be used for extra compression, and `--threads=...` option can be used to specify additional threads for faster compression.

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
${JBROWSE_SOURCE_DIR}/bin/generate-names.pl --tracks genes,variants --hashBits 4 --workdir .
```

## Running JBrowse on Odyssey using Open OnDemand

A JBrowse instance can be launched on Odyssey using the Odyssey Open OnDemand web portal ([https://vdi.rc.fas.harvard.edu/](https://vdi.rc.fas.harvard.edu/)).
From the menu bar, select *Interactive Apps > JBrowse*.
The web form will be used to submit a SLURM job that runs on an Odyssey compute node (see [https://www.rc.fas.harvard.edu/resources/running-jobs/#Slurm_partitions](https://www.rc.fas.harvard.edu/resources/running-jobs/#Slurm_partitions) for a list of some the commonly-available partitions; lab-specific partitions may also be specified, as well as a comma-separated list of partitions that SLURM may consider when scheduling the job).
The custom JBrowse Singularity container image runs the [lighttpd](https://www.lighttpd.net/) web server, which requires only a single processor core and minimal memory (256 MB is requested for the purpose SLURM job scheduling).

In the "path of a JBrowse data directory" textbox, enter the absolute path to the desired JBrowse data/ directory created previously, then click "Launch".
Once the JBrowse job has been scheduled and the JBrowse instance started on a compute node, a "Connect to JBrowse" button will appear; click to open JBrowse in a new tab.

## Support

For assistance with launching JBrowse on Odyssey, please submit a [help ticket](https://portal.rc.fas.harvard.edu/rcrt/submit_ticket).

JBrowse bug reports/feature requests may be submitted to the [JBrowse issue tracker](https://github.com/GMOD/jbrowse/issues) on GitHub.

For general JBrowse questions, see the [JBrowse Gitter chatroom](https://gitter.im/GMOD/jbrowse) and the [gmod-ajax mailing list](https://sourceforge.net/p/gmod/mailman/gmod-ajax/).
