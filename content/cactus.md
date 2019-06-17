Title: Cactus on the FASRC Cluster
Date: 2019-09-19
Author: Nathan Weeks
Category: Software
Tags: Multiple Sequence Alignment, Singularity
Summary: How to run Cactus on the FASRC Cluster

## Introduction

[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) performs Multiple Sequence Alignment on whole genomes.
This guide describes best practices for running Cactus on the FASRC Cluster.

*NOTE: these instructions describe a new protocol for running cactus on the FASRC cluster.
 It has not been tested on large data sets.
 Please [contact us](pages/about) with any questions or feedback.*

## Cactus on the FASRC Cluster

For reliable execution, it is recommended to run Cactus entirely within a provided [Singularity](https://www.rc.fas.harvard.edu/resources/documentation/software/singularity-on-odyssey/) image (`cactus --binariesMode local`) on a single node with enough memory for an MSA of the target genomes.

## Example

The Cactus source distribution includes an example data set that can be retrieved thus:

```
$ curl -LO https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt
```

### Example job script

This job script must be submitted (via sbatch) from a directory on a file system that is accessible from all FASRC cluster compute nodes (e.g., directories prefixed with /n/, such as /n/scratchlfs) that has enough free space to accommodate the Cactus jobStore (at least 100 GB recommended) in the event Cactus does not complete execution before the wall time limit is reached.

For larger genomes, modify the sbatch `--time` option to specify a longer wall clock time limit, and `--partition` option to use larger-memory nodes.
As Cactus memory requirements scale approximately quadratically with genome size, and linearly with the number of genomes, a node from the "general" or "bigmem" partitions (or similar) may be required.

Cactus will automatically limit the number of tasks based on the number of processor cores on the compute node.

    :::sh
    #!/bin/sh
    # allow use of all the memory on the node (Cactus will automatically limit
    # concurrent tasks based on the number of processor cores)
    #SBATCH --nodes=1
    #SBATCH --mem=0
    #SBATCH --exclusive
    # Customize --time --partition as appropriate (general has larger memory & direct-attached storage)
    #SBATCH --time=0:60:00
    #SBATCH --partition=test
    # Terminate Cactus job step (srun command) ~10 minutes (600 seconds) before the --time limit
    # to allow the jobStore to be saved
    #SBATCH --signal=TERM@600
    
    set -o nounset -o errexit -o xtrace
    
    ########################################
    # parameters
    ########################################
    readonly CACTUS_IMAGE=/n/scratchssdlfs/singularity_images/informatics/cactus_2019.09.03-623cfc5.sif 
    readonly SEQFILE=evolverMammals.txt
    readonly OUTPUTHAL=evolverMammals.hal
    # extra options to Cactus
    readonly CACTUS_OPTIONS='--root mr'
    
    # path to squashfs image file of jobstore from previous (incomplete) run, or empty if no previous run
    readonly RESTART_JOBSTORE_IMAGE='' #cactus-jobStore-SLURM_JOB_ID.squashfs
    
    ########################################
    # ... don't modify below here ...
    
    readonly CACTUS_SCRATCH=/scratch/cactus-${SLURM_JOB_ID}

    # Creating
    # Home directory 
    # the toil workDir must be on the same file system as the cactus jobStore
    mkdir -m 700 -p ${CACTUS_SCRATCH}/workDir ${CACTUS_SCRATCH}/tmp ${CACTUS_SCRATCH}/home
    
    if [ "${RESTART_JOBSTORE_IMAGE}" ]
    then
      /usr/bin/time -v /usr/sbin/unsquashfs -d ${CACTUS_SCRATCH}/jobStore "${RESTART_JOBSTORE_IMAGE}"
    fi

    if ! srun -n 1 singularity exec --cleanenv \
                                    --bind ${CACTUS_SCRATCH}:/cactus-scratch \
                                    --bind ${CACTUS_SCRATCH}/tmp:/tmp \
                                    --home ${CACTUS_SCRATCH}/home:/cactus-home ${CACTUS_IMAGE} \
         cactus ${CACTUS_OPTIONS} ${RESTART_JOBSTORE_IMAGE:+--restart} --workDir=/cactus-scratch/workDir --binariesMode local /cactus-scratch/jobStore "${SEQFILE}" "${OUTPUTHAL}"
    then
      /usr/bin/time -v /usr/sbin/mksquashfs ${CACTUS_SCRATCH}/jobStore cactus-jobStore-${SLURM_JOB_ID}.squashfs
    fi 

    # The jobStore directory would eventually be purged, but just in case the
    # next job to run on this node needs lots of /scratch...

    rm -rf ${CACTUS_SCRATCH}

## Restarting an incomplete Cactus run

Using the above job script, the Cactus jobStore will be created on the node-local /scratch file system.
In the event Cactus does not finish within the wall time limit specified in the job script (sbatch `--time` option), a [SquashFS](https://wiki.gentoo.org/wiki/SquashFS) file system image file of the toil jobStore will be staged-out to this directory in a file named `cactus-jobStore-${SLURM_JOB_ID}.squashfs`, where `${SLURM_JOB_ID}` is the numeric ID assigned to the job by SLURM.

The Cactus jobStore contains the state of the Cactus workflow.
When Cactus is run with the `--restart` option, it will resume execution of a workflow, skipping finished tasks.
In the above job script, this can be accomplished by setting the `RESTART_JOBSTORE_IMAGE` variable to the pathname of a Cactus jobStore SquashFS file created by a previous execution of the Cactus job script (e.g., `RESTART_JOBSTORE_IMAGE=cactus-jobStore-12345678.squashfs`) before resubmitting the job script with `sbatch`.

---

*Technical comment*: The description of the exit status of the `srun` command implies that after SIGTERM has been delivered to the job step, (courtesy of the `sbatch --signal=TERM` directive), the `srun` command will block until the tasks have exited, after which it should be safe to copy the jobStore.
From the `srun` man page:

```
RETURN VALUE
    srun will return the highest exit code of all tasks run or the
    highest signal (with the high-order bit set in an 8-bit integer
    -- e.g. 128 + signal)
```

A single SquashFS file is created instead of recursively copying the contents of the jobStore directory as the jobStore may comprise many files, and writing to a single large file on a parallel file system such as Lustre is generally faster than creating/writing to many smaller files.

The `sbatch --signal=TERM@600` in the above script indicates that the TERM signal should be sent approximately 600 seconds before the end of the specified walltime limit (`sbatch --time=00:60:00`).
This amount may be adjusted based on the anticipated size of the Cactus jobStore.
In one test using a node from the (old) FASRC cluster "test" partition (Odyssey nodes, where each node had 2 x 16-core Intel Xeon(R) CPU E5-2683 v4 @ 2.10GHz), the above `mksquashfs` command compressed a Cactus jobStore at a rate of approximately 50 MB/s, with a compression ratio of approximately 6.5X.

---

## About the Cactus Singularity Image

The Cactus Singularity image was generated from the following Dockerfile (specifying the root of a Cactus git working tree as the build context) using [docker2singularity](https://github.com/singularityware/docker2singularity):

```
FROM ubuntu:16.04 AS builder

RUN apt-get update && apt-get install --no-install-recommends -y git gcc g++ build-essential python-dev zlib1g-dev libkyototycoon-dev libtokyocabinet-dev libkyotocabinet-dev libbz2-dev pkg-config python-pip python-setuptools python-wheel

ENV kyotoTycoonIncl -I/usr/include -DHAVE_KYOTO_TYCOON=1
ENV kyotoTycoonLib -L/usr/lib -Wl,-rpath,/usr/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++
RUN mkdir -p /home/cactus

COPY . /home/cactus

RUN cd /home/cactus && make clean && rm -f /home/cactus/submodules/hdf5/bin/h5c++
RUN cd /home/cactus && make

RUN pip install --upgrade --prefix=/home/cactus/.local toil==3.20.0
# https://github.com/ComparativeGenomicsToolkit/cactus/pull/79
RUN cd /home/cactus && sed -i.bak 's/networkx>=2,<3/networkx==2.2/' setup.py && pip install --upgrade --prefix=/home/cactus/.local .

FROM ubuntu:16.04
RUN apt-get update && apt-get install --no-install-recommends -y kyototycoon libtokyocabinet9 python zlib1g bzip2 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /home/cactus/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/sonLib/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/cactus2hal/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/.local/bin/* /usr/local/bin/
# /usr/local/lib/python2.7/site-packages not in sys.path
COPY --from=builder /home/cactus/.local/lib/python2.7/site-packages/ /usr/local/lib/python2.7/dist-packages/
```
