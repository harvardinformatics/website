Title: MAKER on the FASRC Cluster
Date: 2019-06-18
Modified: 2021-01-08
Author: Nathan Weeks
Category: Software
Tags: Genome Annotation, MAKER
Summary: How to run MAKER on the FASRC Cluster

## Introduction

[MAKER](http://www.yandell-lab.org/software/maker.html)[^maker2]<sup>,</sup>[^maker] is a popular genome annotation pipeline for both prokaryotic and eukaryotic genomes.
This guide describes best practices for running MAKER on the FASRC cluster to maximize performance.
For additional MAKER background and examples, see [this tutorial](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018).
Most tutorial examples can be run on an compute node in an [interactive job](https://docs.rc.fas.harvard.edu/kb/running-jobs/#Interactive_jobs_and_srun), prefixing any MAKER commands with `singularity exec --cleanenv ${MAKER_IMAGE}`.
For general MAKER support, see the [maker-devel Google Group](https://groups.google.com/forum/#!forum/maker-devel).

## MAKER on the FASRC Cluster

MAKER is run on the FASRC cluster using a provided [Singularity](https://docs.rc.fas.harvard.edu/kb/singularity-on-the-cluster/) image.
This image was created from the MAKER [Biocontainers](https://biocontainers.pro) image (which was automatically generated from the corresponding [Bioconda package](https://bioconda.github.io/recipes/maker/README.html)), and bundles [RepBase RepeatMasker edition](https://www.girinst.org/repbase/update/) in addition to the default [Dfam 3.2](https://xfam.wordpress.com/2020/07/09/dfam-3-2-release/)[^dfam] library.

## Prerequisites

1. Create the empty MAKER [control files](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained) by running the following [interactive job](https://docs.rc.fas.harvard.edu/kb/running-jobs/#Interactive_jobs_and_srun) from a FAS RC login node (as Singularity is not installed on the FAS RC login nodes):

        :::sh
        srun -p test,serial_requeue,shared sh -c 'singularity exec --cleanenv /n/singularity_images/informatics/maker/maker:2.31.11-repbase.sif maker -CTL'

    This results in 3 files:

    * **maker_opts.ctl** (*required*: modify this file)
    * **maker_exe.ctl** (*do not* modify this file)
    * **maker_bopts.ctl** (*optionally* modify this file)

2. In **maker_opts.ctl**:
    * *(Required)* If not using RepeatMasker, change:

          `model_org=all`
          
          to
          
          `model_org= `
       
          If using RepeatMasker, change `model_org=all` to an appropriate family/genus/species (or other taxonomic rank).
          The [famdb.py](https://github.com/Dfam-consortium/FamDB#famdbpy) utility can be used to query the combined Dfam/Repbase repeat library:


            :::sh
            $ srun --pty test,serial_requeue,shared singularity shell --cleanenv /n/singularity_images/informatics/maker/maker:2.31.11-repbase.sif
            ...
            Singularity> /usr/local/share/RepeatMasker/famdb.py -i /usr/local/share/RepeatMasker/Libraries/RepeatMaskerLib.h5 names Heliconius
            Exact Matches
            =============
            33416 Heliconius (scientific name)

            Non-exact Matches
            =================
            33418 Heliconius ethilla (scientific name), Heliconius ethilla (Godart, 1819) (authority)
            33419 Heliconius numata (scientific name), Heliconius numata (Cramer, 1780) (authority)
            ...
            Singularity> /usr/local/share/RepeatMasker/famdb.py -i /usr/local/share/RepeatMasker/Libraries//RepeatMaskerLib.h5  lineage --ancestors --descendants 'Heliconius melpomene'
            ...
            ...
            └─33416 Heliconius [538]
              └─34740 Heliconius melpomene [75]
                ├─171916 Heliconius melpomene rosina [0]
                ├─171917 Heliconius melpomene melpomene [88]
            ...
            Singularity> exit

          In this case, to use the species "Heliconius melpomene", specify `model_org=heliconius_melpomene` (i.e., case insensitive, replace spaces with underscores)

    * *(Recommended)* Change `max_dna_len=100000` to `max_dna_len=300000` to increase the length of the segments that the reference sequence is partitioned into for sequence alignment.
      This reduces the number of files created during MAKER execution, lessening the file metadata load on the [FASRC scratch file system](https://docs.rc.fas.harvard.edu/kb/policy-scratch/), which is one of the main constraints for MAKER scalability.
    * Other options may need to be adjusted (e.g., `split_hit=10000` corresponds to the longest expected intron length).
      See [this table](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4286374/table/T1/)[^maker] for a description of additional relevant MAKER options in maker_opts.ctl and maker_bopts.ctl.

## Example Job Script

There are two approaches to running MAKER on the FASRC cluster: (1) entirely with in a container on a single compute node (more reliable, but slower) and  (2) on multiple compute nodes (launched using MPI from outside of the container; susceptible to conflicts with the user environment).
Either example job script must be submitted (via sbatch) from a directory on a file system that is mounted on all compute nodes (e.g., directories prefixed with /n/, such as /n/scratchlfs).
The MAKER datastore directory will be created in the directory this job script is submitted from (named using the reference sequence file name prefix, and ending in *-output).

### Example Single-Compute-Node MAKER Job Script

The single-node approach is considered more robust (though less scalable), and is recommended if your if you have environment variables set in bash startup files (e.g., `LD_LIBRARY_PATH` or Perl-related environment variables) that may interfere with the operation of software in the MAKER container.

    :::sh
    #!/bin/sh
    # Customize --time and --partition as appropriate.
    # --exclusive --mem=0 allocates all CPUs and memory on the node.
    #SBATCH --partition=shared
    #SBATCH --nodes=1
    #SBATCH --mem=0
    #SBATCH --exclusive
    #SBATCH --time=0:30:00

    MAKER_IMAGE=/n/singularity_images/informatics/maker/maker:2.31.11-repbase.sif

    # Submit this job script from the directory with the MAKER control files

    # RepeatMasker setup (if not using RepeatMasker, optionally comment-out these three lines)
    export SINGULARITYENV_LIBDIR=${PWD}/LIBDIR
    mkdir -p LIBDIR
    singularity exec ${MAKER_IMAGE} sh -c 'ln -sf /usr/local/share/RepeatMasker/Libraries/* LIBDIR'

    # singularity options:
    # * --cleanenv : don't pass environment variables to container (except those specified in --env option-arguments)
    # * --no-home : don't mount home directory (if not current working directory) to avoid any application/language startup files
    # Add any MAKER options after the "maker" command
    # * -nodatastore is suggested for Lustre, as it reduces the number of directories created
    # * -fix_nucleotides needed for hsap_contig.fasta example data
    singularity exec --no-home --cleanenv ${MAKER_IMAGE} mpiexec -n ${SLURM_CPUS_ON_NODE} maker -fix_nucleotides -nodatastore

### Example Multi-Compute-Node MAKER Job Script 

In the following job script, MAKER can scale across multiple nodes in the FAS RC cluster by increasing the sbatch `--ntasks` value (which indicates the total number of processor cores, or "CPUs", to allocate across any number of compute nodes).
Increase `--ntasks` may cause the job to take longer to schedule and start.
See FAS RC [Slurm Partitions](https://docs.rc.fas.harvard.edu/kb/running-jobs/#Slurm_partitions) for a description of limits on jobs submitted to available Slurm partitions.

    :::sh
    #!/bin/sh
    # Customize --time, --ntasks, and --partition as appropriate
    #SBATCH --time=0:30:00
    #SBATCH --ntasks=8
    #SBATCH --mem-per-cpu=4g
    #SBATCH --partition=shared

    MAKER_IMAGE=/n/singularity_images/informatics/maker/maker:2.31.11-repbase.sif

    # Submit this job script from the directory with the MAKER control files

    # Remove any environment modules
    module purge

    # Use Intel MPI for the "mpiexec" command
    module load intel/19.0.5-fasrc01 impi/2019.8.254-fasrc01

    # RepeatMasker setup (if not using RepeatMasker, optionally comment-out these three lines)
    mkdir -p LIBDIR
    singularity exec ${MAKER_IMAGE} sh -c 'ln -sf /usr/local/share/RepeatMasker/Libraries/* LIBDIR'
    export LIBDIR=$PWD/LIBDIR

    # Add any MAKER options
    # * the -mpi option is needed to use the host MPI for MAKER in a Singularity container
    # * -nodatastore is suggested for Lustre, as it reduces the number of directories created
    # * -fix_nucleotides needed for hsap_contig.fasta example data
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
However, this script executes very slowly on a parallel file system (e.g., FASRC scratchlfs and holylfs file systems), and the resulting JBrowse data directory is completely unsuitable for visualization when located on a parallel file system due to a large number of small files created.
An in-house customization of this script (`ifxmaker2jbrowse`) has been developed and tuned for parallel file systems.
Execution time of `ifxmaker2jbrowse` is well over an order of magnitude faster than `maker2jbrowse`, and the resulting JBrowse data directory contains tens of files in standard formats usable by other tools (e.g., [bgzip](https://www.htslib.org/doc/bgzip.html)-compressed & [tabix](https://www.htslib.org/doc/tabix.html)-indexed [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), and bgzip-compressed and [samtools](https://www.htslib.org/doc/samtools.html)-faidx-indexed FASTA) instead of tens/hundreds of thousands of JBrowse-specific files.

### ifxmaker2jbrowse

A Singularity image containing the [ifxmaker2jbrowse script](images/ifxmaker2jbrowse) and and all dependencies provided.
The following example job script (submitted from the MAKER datastore directory) demonstrates its use:

    :::sh
    #!/bin/sh
    # Customize --time and --partition as appropriate
    #SBATCH --time=0:60:00
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=8
    #SBATCH --mem=32G
    #SBATCH --partition=shared
    
    # set to the pathname of reference FASTA file specified for the maker_opts.ctl "genome=" option.
    REFERENCE_FASTA=../MY_REF.fa 
    
    # options to bgzip and sort (assuming GNU coreutils sort); these are used for optimizing performance and disk usage
    export SINGULARITYENV_MAKER2JBROWSE_BGZIP_OPTIONS="--threads=${SLURM_CPUS_PER_TASK}"
    export SINGULARITYENV_MAKER2JBROWSE_SORT_OPTIONS="--parallel=${SLURM_CPUS_PER_TASK} --stable --buffer-size=1G --compress-program=gzip"
    
    # if you would like to omit the creation of a compressed/indexed reference FASTA file, and store just the
    # reference sequence the lengths for use in JBrowse, add the `--noseq` option to the following command:
    singularity run --cleanenv /n/singularity_images/informatics/maker/ifxmaker2jbrowse:20210108.sif --bgzip_fasta=${REFERENCE_FASTA} --no_names_index --ds_index *_master_datastore_index.log

The recommended ifxmaker2jbrowse `--no_names_index` option disables the creation of a searchable index of all feature names in JBrowse.
If name-based indexing is desired for select tracks, this can subsequently be done more efficiently (resulting in fewer files) using the following options to the JBrowse [generate-names.pl](https://jbrowse.org/docs/generate-names.pl.html) script (e.g., for the protein2genome and est2genome tracks):

```
# execute from the JBrowse data/ directory
singularity exec --cleanenv /n/singularity_images/informatics/maker/ifxmaker2jbrowse:20210108.sif generate-names.pl --tracks protein2genome,est2genome --hashBits 4 --compress --out .
```

### Running JBrowse on the FASRC Cluster using Open OnDemand

A JBrowse instance can be launched on the FASRC cluster using Open OnDemand instance ([https://vdi.rc.fas.harvard.edu/]()).
From the menu, select Interactive Apps > JBrowse.
In the "path of a JBrowse data directory" textbox, enter the absolute path to the JBrowse data/ directory that was created by the ifxmaker2jbrowse script (in the MAKER datastore directory), then click "Launch".

For more details, see the [JBrowse on the FASRC Cluster]({filename}/jbrowse.md) guide.

## Using the MAKER Singularity Image on Other HPC Clusters

The MAKER Singularity image file was obtained from the [Galaxy Project](https://galaxyproject.org/), which maintains a [large repository of Singularity images](https://depot.galaxyproject.org/singularity/) for use by the [Galaxy platform](https://usegalaxy.org/).
It can be downloaded for use in other HPC environments that support Singularity:

```
$ curl -O /n/singularity_images/informatics/maker/maker:2.31.11-repbase
```

## References

[^maker2]: Holt C, Yandell M. *MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects.* BMC Bioinformatics. 2011 Dec 22;12:491. doi: [10.1186/1471-2105-12-491](https://dx.doi.org/10.1186%2F1471-2105-12-491). PMID: 22192575; PMCID: [PMC3280279](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3280279/).
[^maker]: Campbell MS, Holt C, Moore B, Yandell M. *Genome Annotation and Curation Using MAKER and MAKER-P.* Curr Protoc Bioinformatics. 2014 Dec 12;48:4.11.1-39. doi: [10.1002/0471250953.bi0411s48](https://dx.doi.org/10.1002%2F0471250953.bi0411s48). PMID: 25501943; PMCID: [PMC4286374](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4286374/).
[^dfam]: Hubley R, Finn RD, Clements J, Eddy SR, Jones TA, Bao W, Smit AF, Wheeler TJ. *The Dfam database of repetitive DNA families.* Nucleic Acids Res. 2016 Jan 4;44(D1):D81-9. doi: [10.1093/nar/gkv1272](https://dx.doi.org/10.1093%2Fnar%2Fgkv1272). Epub 2015 Nov 26. PMID: 26612867; PMCID: [PMC4702899](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4702899/).
