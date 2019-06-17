Title: Cactus on Odyssey
Date: 2019-06-17
Author: Nathan Weeks
Category: Tutorials
Tags: Odyssey, Multiple Sequence Alignment
Summary: How to run Cactus on Odyssey

## Introduction

[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) performs Multiple Sequence Alignment on whole genomes.
This guide describes best practices for running Cactus on Odyssey.

## Cactus on Odyssey

For reliable execution, it is recommended to run Cactus entirely within a provided [Singularity](https://www.rc.fas.harvard.edu/resources/documentation/software/singularity-on-odyssey/) image (`cactus --binariesMode local`) on a single node with enough memory for an MSA of the target genomes.
As Cactus memory requirements scale approximately quadratically with genome size, and linearly with the number of genomes, a node from the "bigmem" partition (or similar) may be required.

The provided Cactus Singularity image was created from the Cactus [Biocontainers](https://biocontainers.pro) image, which in turn was generated from the corresponding [Bioconda package](https://bioconda.github.io/recipes/cactus/README.html).

## Example

The Cactus source distribution includes an example data set ("evolverMammals").

Due to a known issue affecting DNS in the Cactus Singularity image, input files specified in the [evolverMammals.txt](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt) seqFile must be downloaded manually:

```
$ mkdir -p evolverMammals
$ curl http://s3-us-west-2.amazonaws.com/jcarmstr-misc/testRegions/evolverMammals/simCow.chr6 > evolverMammals/simCow.chr6
$ curl http://s3-us-west-2.amazonaws.com/jcarmstr-misc/testRegions/evolverMammals/simDog.chr6 > evolverMammals/simDog.chr6
$ curl http://s3-us-west-2.amazonaws.com/jcarmstr-misc/testRegions/evolverMammals/simHuman.chr6 > evolverMammals/simHuman.chr6
$ curl http://s3-us-west-2.amazonaws.com/jcarmstr-misc/testRegions/evolverMammals/simMouse.chr6 > evolverMammals/simMouse.chr6
$ curl http://s3-us-west-2.amazonaws.com/jcarmstr-misc/testRegions/evolverMammals/simRat.chr6 > evolverMammals/simRat.chr6
```

Create the input seqFile (evolverMammals.txt) to have the following contents:

```
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974):0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303):0.032898);

simCow_chr6 evolverMammals/simCow.chr6
simDog_chr6 evolverMammals/simDog.chr6
simHuman_chr6 evolverMammals/simHuman.chr6
simMouse_chr6 evolverMammals/simMouse.chr6
simRat_chr6 evolverMammals/simRat.chr6
```

### Example job script

This job script must be submitted (via sbatch) from a directory on a file system that is accessible from all Odyssey compute nodes (e.g., directories prefixed with /n/, such as /n/scratchlfs) that has enough free space to accommodate the Cactus jobStore (at least 100 GB recommended) in the event Cactus does not complete execution before the wall time limit is reached.

For larger genomes, modify the sbatch `--time` option to specify a longer wall clock time limit, and `--partition` option to use larger-memory nodes.

Cactus will automatically limit the number of tasks based on the number of processor cores on the compute node.

    :::sh
    #!/bin/sh
    # allow use of all the memory on the node (Cactus will automatically limit
    # concurrent tasks based on the number of processor cores)
    #SBATCH --nodes=1
    #SBATCH --mem=0
    #SBATCH --exclusive
    # Customize --time --partition as appropriate
    #SBATCH --time=0:30:00
    #SBATCH --partition=test,shared
    # Terminate Cactus job step (srun command) ~10 minutes (600 seconds) before the --time limit
    # to allow the jobStore to be saved
    #SBATCH --signal=KILL@600
    
    set -o nounset -o errexit -o xtrace
    
    ########################################
    # parameters
    ########################################
    readonly CACTUS_IMAGE=/n/scratchssdlfs/singularity_images/informatics/cactus:2019.03.01--py27hdbcaa40_1.sif 
    readonly SEQFILE=evolverMammals.txt
    readonly OUTPUTHAL=evolverMammals.hal
    # extra options to Cactus
    readonly CACTUS_OPTIONS='--root mr'
    
    # path to tar file of jobstore from previous (incomplete) run, or empty if no previous run
    readonly RESTART_JOBSTORE_TAR_FILE='' # cactus-jobStore-${SLURM_JOB_ID}.tar.gz
    
    ########################################
    # ... don't modify below here ...
    
    readonly JOBSTORE=/scratch/cactus-jobStore-${SLURM_JOB_ID}
    
    if [ "${RESTART_JOBSTORE_TAR_FILE}" ]
    then
      mkdir ${JOBSTORE}
      tar -C ${JOBSTORE} -xzf "${RESTART_JOBSTORE_TAR_FILE}"
    fi
    
    if ! srun -n 1 singularity exec ${CACTUS_IMAGE} cactus ${RESTART_JOBSTORE_TAR_FILE:+--restart} --binariesMode local ${CACTUS_OPTIONS} ${JOBSTORE} "${SEQFILE}" "${OUTPUTHAL}"
    then
      tar -C ${JOBSTORE} -cf - . | gzip -1 > ${JOBSTORE##*/}.tar.gz
    fi 
    
    # The jobStore directory would eventually be purged, but just in case the
    # next job to run on this node needs lots of /scratch...
    rm -rf ${JOBSTORE}

## Restarting an incomplete Cactus run

Using the above job script, the Cactus jobStore will be created on the node-local /scratch file system.
In the event Cactus does not finish within the wall time limit specified in the job script (sbatch `--time` option), a compressed [tar](https://en.wikipedia.org/wiki/Tar_(computing)) file of the toil jobStore will be staged-out to this directory in a file named `cactus-jobStore-${SLURM_JOB_ID}.tar.gz`, where `${SLURM_JOB_ID}` is the numeric ID assigned to the job by SLURM.

The Cactus jobStore contains the state of the Cactus workflow.
When Cactus is run with the `--restart` option, it will resume execution of a workflow, skipping finished tasks.
In the above job script, this can be accomplished by setting the `RESTART_JOBSTORE_TAR_FILE` variable to the pathname of a Cactus jobStore tar file created by a previous execution of the Cactus job script (e.g., `RESTART_JOBSTORE_TAR_FILE=cactus-jobStore-12345678.tar.gz`) before resubmitting the job script with `sbatch`.

---

*Technical comment*: The description of the exit status of the `srun` command implies that after SIGKILL has been delivered to the job steps (courtesy of the `sbatch --signal=KILL` directive), the `srun` command will block until the tasks have exited, after which it should be safe to tar the jobStore.
From the `srun` man page:

```
RETURN VALUE
    srun will return the highest exit code of all tasks run or the
    highest signal (with the high-order bit set in an 8-bit integer
    -- e.g. 128 + signal)
```

A single tar file is created instead of recursively copying the contents of the jobStore directory as the jobStore may comprise many files, and writing to a single large file on a parallel file system such as Lustre is generally faster than creating/writing to many smaller files.

---

## Using the Cactus Singularity image on other HPC clusters

The Cactus Singularity image file can be transferred to another HPC cluster, or the Docker image can be downloaded from the [quay.io](https://quay.io) container registry and converted into a Singularity image using the `singularity pull` command on the target system:

```
singularity pull cactus:2019.03.01--py27hdbcaa40_1.sif docker://quay.io/biocontainers/cactus:2019.03.01--py27hdbcaa40_1.sif 
```
