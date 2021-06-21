Title: Cactus on the FASRC Cluster
Date: 2020-07-30
Modified: 2021-06-21
Author: Nathan Weeks
Category: Software
Tags: Multiple Sequence Alignment, Singularity
Summary: How to run Cactus on the FASRC Cluster

## Introduction

[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) performs Multiple Sequence Alignment on whole genomes.
This guide describes best practices for running Cactus on the [FASRC Cluster](https://www.rc.fas.harvard.edu/cluster/).

*NOTE: these instructions describe a new protocol for running cactus on the FASRC cluster.
 It has not been tested on large data sets.
 Please [contact us](pages/about) with any questions or feedback.*

## Cactus on the FASRC Cluster

For reliable execution, it is recommended to run Cactus entirely within a provided [Singularity](https://docs.rc.fas.harvard.edu/kb/singularity-on-the-cluster/) image (`cactus --binariesMode local`), as illustrated below, on a single node with enough memory for an MSA of the target genomes.

## Example

The Cactus source distribution includes an example data set that can be retrieved thus:

```
$ curl -LO https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt
```

### Example job script

This job script must be submitted (via [sbatch](https://docs.rc.fas.harvard.edu/kb/running-jobs/#Submitting_batch_jobs_using_the_sbatch_command)) from a directory on a file system that is accessible from all FASRC cluster compute nodes (e.g., directories prefixed with /n/, such as a subdirectory of the scratch directory identified by the [$SCRATCH](https://docs.rc.fas.harvard.edu/kb/policy-scratch/) environment variable) that has enough free space to accommodate the Cactus [job store](https://toil.readthedocs.io/en/latest/running/introduction.html#job-store), which stores the state of the workflow.
Depending on the number, size, and phylogeny of the genomes, the job store may require hundreds of gigabytes of storage.

To minimize the performance impact of running Cactus on a parallel file system, this example job script creates the Cactus job store in an image file (`jobStore.img`) that contains an EXT3 file system used as a [persistent overlay](https://sylabs.io/guides/3.5/user-guide/persistent_overlays.html) for the Singularity container.
The job store image file is created as a [sparse file](https://en.wikipedia.org/wiki/Sparse_file), so it initially consumes far less disk space (as indicated by the `du` command) than its maximum size (2T):

```
$ ls -lh jobStore.img
-rw-rw----+ 1 user group 2.0T Dec  5 16:22 jobStore.img
$ du -h jobStore.img
34G     jobStore.img
```

As Cactus memory requirements scale approximately quadratically with genome size, and linearly with the number of genomes, a node from the "bigmem" partitions (or similar) may be required.
Cactus will automatically limit the number of concurrent tasks based on the number of processor cores detected on the compute node.

    :::sh
    #!/bin/sh
    #SBATCH --nodes=1
    # allow use of all the memory on the node
    #SBATCH --mem=0
    # request all CPU cores on the node
    #SBATCH --exclusive
    # Customize --time --partition as appropriate
    #SBATCH --time=0:30:00
    #SBATCH --partition=test
    
    set -o nounset -o errexit -o xtrace
    
    ########################################
    # parameters
    ########################################
    readonly CACTUS_IMAGE=/n/singularity_images/informatics/cactus/cactus_v2.0.1.sif
    readonly JOBSTORE_IMAGE=jobStore.img # cactus jobStore; will be created if it doesn't exist
    readonly SEQFILE=evolverMammals.txt
    readonly OUTPUTHAL=evolverMammals.hal
    # extra options to Cactus
    readonly CACTUS_OPTIONS='--root mr' # NOTE: specific to evolverMammals.txt; change/remove for other input seqFile
    
    ########################################
    # ... don't modify below here ...
    
    readonly CACTUS_SCRATCH=/scratch/cactus-${SLURM_JOB_ID}
    
    if [ ! -e "${JOBSTORE_IMAGE}" ]
    then
      restart=''
      mkdir -p -m 777 ${CACTUS_SCRATCH}/upper ${CACTUS_SCRATCH}/work
      truncate -s 2T "${JOBSTORE_IMAGE}"
      singularity exec ${CACTUS_IMAGE} mkfs.ext3 -d ${CACTUS_SCRATCH} "${JOBSTORE_IMAGE}"
    else
      restart='--restart'
    fi
    
    # Use empty /tmp directory in the container
    mkdir -m 700 -p ${CACTUS_SCRATCH}/tmp
    
    # the toil workDir must be on the same file system as the cactus jobStore
    singularity exec --overlay ${JOBSTORE_IMAGE} ${CACTUS_IMAGE} mkdir -p /cactus/workDir
    srun -n 1 /usr/bin/time -v singularity exec --cleanenv \
                               --overlay ${JOBSTORE_IMAGE} \
                               --bind ${CACTUS_SCRATCH}/tmp:/tmp \
                               --env PYTHONNOUSERSITE=1 \
                               ${CACTUS_IMAGE} \
      cactus ${CACTUS_OPTIONS-} ${restart-} --workDir=/cactus/workDir --binariesMode local /cactus/jobStore "${SEQFILE}" "${OUTPUTHAL}"
    
    # /tmp would eventually be purged, but just in case the
    # next job to run on this node needs lots of /space...
    
    rm -rf ${CACTUS_SCRATCH} jobStore.img

## Restarting an incomplete Cactus run

In the event Cactus does not finish within the wall time limit specified in the job script (sbatch `--time` option) and is terminated by SLURM, simply resubmit the job script with the `sbatch` command.

## Troubleshooting

**Problem:** If spaces appear in the FASTA defline, an error like the following will occur:

> RuntimeError: The fasta header 'CM000102.5 Gallus gallus chromosome 10' contains spaces or tabs. These characters will cause issues in space-separated formats like MAF, and may not function properly when viewed in a browser. Please remove these characters from the input headers and try again.

**Solution:** The following sed command can be used to remove the first space (and all subsequent characters) from FASTA deflines (e.g., in the above example, `>CM000102.5 Gallus gallus chromosome 10` becomes `CM000102.5`):

```
sed '/^>/s/ .*//'  my.fa > my-mod-defline.fa
```

## Working with HAL files

Cactus outputs multiple sequence alignment and ancestral reconstruction information in a [Hierarchical Alignment (HAL)](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md) file.
[HAL Tools](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md#hal-tools) is a suite of command-line utilities for analyzing HAL files and converting to/from different formats.

HAL Tools is installed in the Cactus Singularity image.
For example, the following commands run [halValidate](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md#halvalidate) and [halStats](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md#halstats) on the example output evolverMammals.hal:

*NOTE: on the FAS RC Cluster, Singularity is not available on login nodes, and these commands must be run [within a batch or interactive environment](https://docs.rc.fas.harvard.edu/kb/singularity-on-the-cluster/#Singularity_on_the_cluster) on a compute node*

    :::sh
    singularity exec --cleanenv /n/singularity_images/informatics/cactus/cactus_v2.0.1.sif halValidate evolverMammals.hal
    singularity exec --cleanenv /n/singularity_images/informatics/cactus/cactus_v2.0.1.sif halStats evolverMammals.hal

## About the Cactus Singularity Image

The Cactus Singularity image was generated from the [cactus:v2.0.1 release](https://github.com/ComparativeGenomicsToolkit/cactus/releases/tag/v2.0.1) using the `singularity pull` command:

```
singularity pull --disable-cache docker://quay.io/comparative-genomics-toolkit/cactus:v2.0.1
```
