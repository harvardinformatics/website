Title: MAKER on Odyssey
Date: 2019-06-18
Author: Nathan Weeks
Category: Software
Tags: Odyssey, Genome Annotation, MAKER
Summary: How to run MAKER on Odyssey

## Introduction

[MAKER](http://www.yandell-lab.org/software/maker.html) is a popular genome annotation pipeline for both prokaryotic and eukaryotic genomes.
This guide describes best practices for running MAKER on Odyssey to maximize performance.
For additional MAKER background and examples, see [this tutorial](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018).
Most tutorial examples can be run on an Odyssey compute node in an [interactive job](https://www.rc.fas.harvard.edu/resources/running-jobs/#Interactive_jobs_and_srun), prefixing any MAKER commands with `singularity exec ${MAKER_IMAGE}`.
For general MAKER support, see the [maker-devel Google Group](https://groups.google.com/forum/#!forum/maker-devel).

## MAKER on Odyssey

MAKER is run on Odyssey using a provided [Singularity](https://www.rc.fas.harvard.edu/resources/documentation/software/singularity-on-odyssey/) image.
This image was created from the MAKER [Biocontainers](https://biocontainers.pro) image (which in turn was generated from the corresponding [Bioconda package](https://bioconda.github.io/recipes/maker/README.html)), and at the time of this writing contains the latest versions of most bioinformatics software dependencies.

## Prerequisites

These instructions assume the shell variable `MAKER_IMAGE=/n/scratchssdlfs/singularity_images/informatics/maker:2.31.10--pl526_14.sif` is set.

1. Create the empty MAKER [control files](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained) via `maker -CTL` (required for maker_exe.ctl to generate pathnames for applications inside the container; existing maker_opts.ctl and maker_bopts.ctl files may be used if desired).
   Singularity is not available on Odyssey login nodes, however, `srun` may be used to run `maker -CTL` in a Singularity container on compute node and generate the MAKER control files in the current working directory:
   ```
   srun singularity exec --cleanenv ${MAKER_IMAGE} maker -CTL
   ```

2. In maker_opts.ctl, set `max_dna_len=999999999` to avoid splitting reference sequence contigs into smaller segments for sequence alignment.
   This lessens the file metadata load on the parallel file system (scratchlfs or holylfs on Odyssey), which is the main constraint for MAKER scalability on Odyssey.
   Note that with a test data set, memory usage per process remained below the `--mem-per-cpu=4g` limit suggested in the SLURM job script below, but YMMV.

3. *(Optional; applicable to eukaryotic genomes only)* Download [RepBase RepeatMasker Edition](https://www.girinst.org/) to use as a supplemental repeat database for [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) (license required == $$).
   See the comments of the sample job script below for where to copy the RepeatMaskerLib.embl file.
   Note that RepeatMasker is bundled with [Dfam 3.0](https://dfam.org/home).

## Example Job Script

This job script must be submitted (via sbatch) from a directory on a file system that is mounted on all Odyssey compute nodes (e.g., directories prefixed with /n/, such as /n/scratchlfs).
The MAKER datastore directory will be created in the directory this job script is submitted from (named using the reference sequence file name prefix, and ending in *-output).

    :::sh
    #!/bin/sh
    # Customize --time, --ntasks, and --partition as appropriate
    #SBATCH --time=0:30:00
    # --ntasks should be <= 100 on Odyssey!!!
    #SBATCH --ntasks=8
    #SBATCH --mem-per-cpu=4g
    #SBATCH --partition=shared
    
    MAKER_IMAGE=/n/scratchssdlfs/singularity_images/informatics/maker:2.31.10--pl526_14.sif
    
    # NOTE: empty MAKER control files can be generated using the command:
    #       singularity exec ${MAKER_IMAGE} maker -CTL 
    # This will be needed at least for the maker_exe.ctl file, which has the paths to executables in the container.
    # Otherwise, existing maker_bopts.ctl and maker_opts.ctl should be usable.
    
    # Submit this job script from the directory with the MAKER control files
    
    # Remove any environment variables
    module purge
    
    # Use Intel MPI for the "mpiexec" command
    module load intel/17.0.4-fasrc01 impi/2017.2.174-fasrc01
    
    # Optional repeat masking (if not using RepeatMasker, comment-out these three lines)
    export REPEATMASKER_LIB_DIR=$PWD/REPEATMASKER_LIB_DIR
    mkdir -p "${REPEATMASKER_LIB_DIR}"
    singularity exec ${MAKER_IMAGE} sh -c "ln -sf /usr/local/share/RepeatMasker/Libraries/* ${REPEATMASKER_LIB_DIR}"
    # If RepBase RepeatMasker Edition has been downloaded, it should be copied into this directory:
    #   cp /path/to/RepeatMaskerLib.embl ${REPEATMASKER_LIB_DIR}
    
    # Add any MAKER options (-fix_nucleotides needed for hsap_contig.fasta example data)
    # * the -mpi option is needed to use MPI with MAKER in a Singularity container
    # * -nodatastore is suggested for Lustre, as it reduces the number of directories created
    mpiexec -n ${SLURM_NTASKS} singularity exec ${MAKER_IMAGE} maker -mpi -fix_nucleotides -nodatastore

---
*NOTE*: MAKER will emit the following warnings during execution; they can be ignored:
```
Possible precedence issue with control flow operator at /usr/local/lib/site_perl/5.26.2/Bio/DB/IndexedBase.pm line 805.

df: Warning: cannot read table of mounted file systems: No such file or directory
```
---

## Visualizing in JBrowse

[JBrowse](https://jbrowse.org/) can be used to visualize MAKER-generated annotation, RNA/protein evidence sequence alignments, and RepeatMasker-masked regions in the context of the reference genome.

JBrowse provides a script `maker2jbrowse` that automatically exports a MAKER datastore to a JBrowse data directory that can be directly visualized in JBrowse.
However, this script executes very slowly on a parallel file system (such as Odyssey's scratchlfs and holylfs), and the resulting JBrowse data directory is completely unsuitable for visualization when located on a parallel file system due to a large number of small files created.
An in-house customization of this script (`maker2jbrowse-odyssey`) has been developed and tuned for parallel file systems.
Execution time of `maker2jbrose-odyssey` is well over an order of magnitude faster than `maker2jbrowse`, and the resulting JBrowse data directory contains tens of files in standard formats usable by other tools (e.g., [bgzip](https://www.htslib.org/doc/bgzip.html)-compressed & [tabix](https://www.htslib.org/doc/tabix.html)-indexed [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), and bgzip-compressed and [samtools](https://www.htslib.org/doc/samtools.html)-faidx-indexed FASTA) instead of tens/hundreds of thousands of JBrowse-specific files.

### maker2jbrowse-odyssey setup

Create a [conda environment](https://informatics.fas.harvard.edu/python-on-odyssey.html) with the [JBrowse Perl libraries and command-line utilities](https://bioconda.github.io/recipes/jbrowse/README.html), [samtools](https://bioconda.github.io/recipes/samtools/README.html), and [htslib](https://bioconda.github.io/recipes/htslib/README.html), then copy the [maker2jbrowse-odyssey script](images/maker2jbrowse-odyssey) into the ${JBROWSE_SOURCE_DIR}/bin directory and ensure it is executable:

```
$ module load Anaconda3/5.0.1-fasrc02
$ conda create -n jbrowse-utils jbrowse samtools htslib
$ source activate jbrowse-utils
(jbrowse-utils) $ cp /path/to/maker2jbrowse-odyssey ${JBROWSE_SOURCE_DIR}/bin
(jbrowse-utils) $ chmod +x ${JBROWSE_SOURCE_DIR}/bin/maker2jbrowse-odyssey
```

### maker2jbrowse-odyssey execution

The following example job script (submitted from the MAKER datastore directory) demonstrates the use of the maker2jbrowse-odyssey script:

    :::sh
    #!/bin/sh
    # Customize --time and --partition as appropriate
    #SBATCH --time=0:60:00
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=8
    #SBATCH --mem=32G
    #SBATCH --partition=shared
    
    module load Anaconda3/5.0.1-fasrc02
    source activate jbrowse-utils
    
    # set to the pathname of reference FASTA file specified for the maker_opts.ctl "genome=" option.
    REFERENCE_FASTA=../MY_REF.fa 
    
    # options to bgzip and sort (assuming GNU coreutils sort); these are used for optimizing performance and disk usage
    export MAKER2JBROWSE_BGZIP_OPTIONS="--threads=${SLURM_CPUS_PER_TASK}"
    export MAKER2JBROWSE_SORT_OPTIONS="--parallel=${SLURM_CPUS_PER_TASK} --stable --buffer-size=1G --compress-program=gzip"
    
    # if you would like to omit the creation of a compressed/indexed reference FASTA file, and store just the
    # reference sequence the lengths for use in JBrowse, add the `--noseq` option to the following command:
    ${JBROWSE_SOURCE_DIR}/bin/maker2jbrowse-odyssey --bgzip_fasta=${REFERENCE_FASTA} --no_names_index --ds_index *_master_datastore_index.log

The recommended maker2jbrowse-odyssey `--no_names_index` option disables the creation of a searchable index of all feature names in JBrowse.
If name-based indexing is desired for select tracks, this can subsequently be done more efficiently (resulting in fewer files) using the following options to the JBrowse [generate-names.pl](https://jbrowse.org/docs/generate-names.pl.html) script (e.g., for the protein2genome and est2genome tracks):

```
# execute from the JBrowse data/ directory
(jbrowse-utils) $ ${JBROWSE_SOURCE_DIR}/bin/generate-names.pl --tracks protein2genome,est2genome --hashBits 4 --compress --out .
```

### Running JBrowse on Odyssey using Open OnDemand

A JBrowse instance can be launched on Odyssey using the Odyssey Open OnDemand instance ([https://vdi.rc.fas.harvard.edu/]()).
From the menu, select Interactive Apps > JBrowse.
In the "path of a JBrowse data directory" textbox, enter the absolute path to the JBrowse data/ directory that was created by the maker2jbrowse-odyssey script (in the MAKER datastore directory), then click "Launch".

For more details, see the [JBrowse on Odyssey guide]({filename}/jbrowse.md).

## Using the MAKER Singularity Image on Other HPC Clusters

The MAKER Singularity image file can be transferred to another HPC cluster, or the Docker image can be downloaded from the [quay.io](https://quay.io) container registry and converted into a Singularity image using the `singularity pull` command on the target system:

```
singularity pull maker:2.31.10--pl526_14.sif docker://quay.io/biocontainers/maker:2.31.10--pl526_14
```
