Title: Cactus on the FASRC Cluster
Date: 2019-09-19
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

For reliable execution, it is recommended to run Cactus entirely within a provided [Singularity](https://www.rc.fas.harvard.edu/resources/documentation/software/singularity-on-odyssey/) image (`cactus --binariesMode local`) on a single node with enough memory for an MSA of the target genomes.

## Example

The Cactus source distribution includes an example data set that can be retrieved thus:

```
$ curl -LO https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt
```

### Example job script

This job script must be submitted (via [sbatch](https://www.rc.fas.harvard.edu/resources/running-jobs/#Submitting_batch_jobs_using_the_sbatch_command)) from a directory on a file system that is accessible from all FASRC cluster compute nodes (e.g., directories prefixed with /n/, such as a subdirectory of the scratch directory identified by the [$SCRATCH](https://www.rc.fas.harvard.edu/policy-scratch/) environment variable) that has enough free space to accommodate the Cactus [job store](https://toil.readthedocs.io/en/latest/running/introduction.html#job-store), which stores the state of the workflow.
Depending on the number, size, and phylogeny of the genomes, the job store may require hundreds of gigabytes of storage.

To minimize the performance impact of running Cactus on a parallel file system, this example job script creates the Cactus job store in an image file (`jobStore.img`) that contains an EXT3 file system used as a [persistent overlay](https://sylabs.io/guides/3.5/user-guide/persistent_overlays.html) for the Singularity container.
The job store image file is created as a [sparse file](https://en.wikipedia.org/wiki/Sparse_file), so it initially consumes far less disk space (as indicated by the `du` command) than its maximum size (1T):

```
$ ls -lh jobStore.img
-rw-rw----+ 1 user group 1.0T Dec  5 16:22 jobStore.img
$ du -h jobStore.img
18G     jobStore.img
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
    readonly CACTUS_IMAGE=/n/singularity_images/informatics/cactus/cactus:2019-11-29.sif
    readonly JOBSTORE_IMAGE=jobStore.img # cactus jobStore; will be created if it doesn't exist
    readonly SEQFILE=evolverMammals.txt
    readonly OUTPUTHAL=evolverMammals.hal
    # extra options to Cactus
    readonly CACTUS_OPTIONS='--root mr'
    
    ########################################
    # ... don't modify below here ...
    
    readonly CACTUS_SCRATCH=/scratch/cactus-${SLURM_JOB_ID}
    
    if [ ! -e "${JOBSTORE_IMAGE}" ]
    then
      restart=''
      mkdir -m 700 -p ${CACTUS_SCRATCH}/upper
      truncate -s 1T "${JOBSTORE_IMAGE}"
      singularity exec /n/singularity_images/informatics/cactus/e2fsprogs:1.45.2.sif mkfs.ext3 -d ${CACTUS_SCRATCH} "${JOBSTORE_IMAGE}"
    else
      restart='--restart'
    fi
    
    # Use empty home & /tmp directories in the container (to avoid, e.g., pip-installed packages in ~/.local)
    mkdir -m 700 -p ${CACTUS_SCRATCH}/home ${CACTUS_SCRATCH}/tmp
    
    # the toil workDir must be on the same file system as the cactus jobStore
    singularity exec --overlay ${JOBSTORE_IMAGE} ${CACTUS_IMAGE} mkdir -p /cactus/workDir
    srun -n 1 singularity exec --cleanenv \
                               --overlay ${JOBSTORE_IMAGE} \
                               --bind ${CACTUS_SCRATCH}/tmp:/tmp \
                               --home ${CACTUS_SCRATCH}/home:/cactus/home ${CACTUS_IMAGE} \
      cactus ${CACTUS_OPTIONS} ${restart-} --workDir=/cactus/workDir --binariesMode local /cactus/jobStore "${SEQFILE}" "${OUTPUTHAL}"
    
    # /tmp would eventually be purged, but just in case the
    # next job to run on this node needs lots of /space...
    
    rm -rf ${CACTUS_SCRATCH}

## Restarting an incomplete Cactus run

In the event Cactus does not finish within the wall time limit specified in the job script (sbatch `--time` option) and is terminated by SLURM, simply resubmit the job script with the `sbatch` command.

## Troubleshooting

**Problem:** If spaces appear in the FASTA defline, an error like the following will occur:

> RuntimeError: The fasta header 'CM000102.5 Gallus gallus chromosome 10' contains spaces or tabs. These characters will cause issues in space-separated formats like MAF, and may not function properly when viewed in a browser. Please remove these characters from the input headers and try again.

**Solution:** The following sed command can be used to remove the first space (and all subsequent characters) from FASTA deflines (e.g., in the above example, `>CM000102.5 Gallus gallus chromosome 10` becomes `CM000102.5`):

```
sed '/^>/s/ .*//'  my.fa > my-mod-defline.fa
```

---

## About the Cactus Singularity Image

The Cactus Singularity image was generated from the following Dockerfile (specifying the root of a Cactus git working tree as the build context) using [docker2singularity](https://github.com/singularityware/docker2singularity):

```
# https://github.com/ComparativeGenomicsToolkit/cactus/issues/94
FROM ubuntu:xenial-20191108 AS builder

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

FROM ubuntu:xenial-20191108
# choose appropriate fast mirror
RUN sed -i'' 's#http://archive.ubuntu.com/ubuntu/#http://mirror.math.princeton.edu/pub/ubuntu/#' /etc/apt/sources.list
RUN apt-get update && apt-get install --no-install-recommends -y kyototycoon libtokyocabinet9 python python-pytest zlib1g bzip2 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /home/cactus/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/sonLib/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/cactus2hal/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/.local/bin/* /usr/local/bin/
# /usr/local/lib/python2.7/site-packages not in sys.path
COPY --from=builder /home/cactus/.local/lib/python2.7/site-packages/ /usr/local/lib/python2.7/dist-packages/

# https://github.com/ComparativeGenomicsToolkit/cactus/issues/102#issuecomment-540781259
RUN sed -i'' -n '/config.writeXML/d; p' /usr/local/lib/python2.7/dist-packages/cactus/progressive/cactus_createMultiCactusProject.py
```
