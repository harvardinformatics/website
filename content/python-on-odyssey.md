Title: Python on Odyssey
Date: 2018-05-23
Author: Aaron Kitzmiller
Category: Software
Tags: Unix, Python, Odyssey
Summary: Using Python on the Odyssey cluster

[TOC]

*Updated for Odyssey 3*

This page describes the use of Python on the Odyssey cluster.  See also the [RC Website](https://www.rc.fas.harvard.edu/resources/documentation/software-on-odyssey/python/) and the [Practical Python on Odyssey tutorial]({filename}/practical-python-on-odyssey.md)

## The Odyssey system Python is 2.7.5 and you can't do much with it
As of this writing, Odyssey uses a CentOS 7 operating system, which includes Python version 2.7.5.  In addition to being somewhat old, the site-packages directory is owned by root.  As a result, any attempt to install a python package will be rebuffed.


## On Odyssey, Python should be used within an Anaconda environment
Python is an extremely popular interpreted programming language in bioinformatics that has an array of excellent packages for scientific computing. It is not uncommon to see module systems on academic clusters provide some of the major Python packages like scipy and numpy. However, in our experience, this leads to complex PYTHONPATHs that can create incompatible environments.

Use of the [Anaconda distribution from Continuum Analytics](https://docs.continuum.io/anaconda/) makes local control of Python package sets much easier. In addition to a large set of pre-installed packages, it is easy to create a local environment in your home directory and use that to define your environment.

You can use Python and Anaconda on Odyssey by running:

    :::bash
    $ module load python/2.7.14-fasrc01

Loading this python module will load the 5.0.1 version of the Anaconda distribution.

Anaconda has a concept of environments that can be used to manage alternative package sets and versions. This analogous to the environments that can be setup with virtualenv. For example, if you want newer versions of some packages in the default environment, you can make a new environment with your customizations.

First, load the base environment module, then create a new environment by cloning the Anaconda root distribution (substitute `ody` with whatever name you wish) and activate it.

{% include 'anaconda-setup.md' %}

If you want to use this environment all the time, add the above line to your ~/.bashrc (or other appropriate shell config file) after the line that loads the module.

You can see where your clone is installed by doing a `which python` when it is activated.


To stop using the custom environment, run:

    :::bash
    $ source deactivate

## You can share a clone with collaborators by installing with a full path (`-p`)
When creating a clone with the `-n` switch, you are "naming" the clone and installing it under your home directory. On Odyssey, this is typically under `$HOME/.conda/envs`.

To share a clone with members of your laboratory or a collaborator, you can install it under a shared location using the full path (`-p`):

    :::bash
    $ module load python/2.7.14-fasrc01
    $ conda create -p /n/pi_lab/software/envs/picrust --clone $PYTHON_HOME
    $ source activate /n/pi_lab/software/envs/picrust
    (/n/pi_lab/software/envs/picrust) $

## You can activate a clone by updating your PATH
The standard way to activate a clone is to use the `source activate` command after loading the appropriate python module.  This only does 2 things: 1) add the bin directory of your clone to your PATH variable and 2) setup the prompt to include the name of the clone.  The ability to `source deactivate` is also handy.

However, if you only use a single clone and don't need the prompt modification, you can simply add the bin directory of your clone to your PATH variable directly.  Then you don't even need to load the python module.

    :::bash
    $ export PATH=$HOME/.conda/envs/ody/bin:$PATH
    $ which python
    ~/.conda/envs/ody/bin/python

No modifications of PYTHONPATH are necessary.

## Within a conda clone, you can install python packages with `conda install`.
In addition to environment management, Anaconda provides package installation through the `conda install` command.  This works just like pip for pure Python packages.  However conda packages can go above and beyond just Python and include compiled C libraries and executables.  For example, the conda install of PySam includes the samtools C libraries:

    :::bash
    (ody) $ conda install pysam
    Fetching package metadata ...............
    Solving package specifications: .

    Package plan for installation in environment /n/home_rc/akitzmiller/.conda/envs/ody:

    The following NEW packages will be INSTALLED:

        bcftools: 1.5-0                bioconda
        bzip2:    1.0.6-1              conda-forge
        htslib:   1.5-0                bioconda
        pysam:    0.11.2.2-htslib1.5_2 bioconda

    The following packages will be UPDATED:

        samtools: 1.4-0                bioconda    --> 1.5-0 bioconda

    Proceed ([y]/n)? y

    bzip2-1.0.6-1. 100% |#############################################################| Time: 0:00:00   2.21 MB/s
    bcftools-1.5-0 100% |#############################################################| Time: 0:00:00   7.60 MB/s
    htslib-1.5-0.t 100% |#############################################################| Time: 0:00:00  14.61 MB/s
    samtools-1.5-0 100% |#############################################################| Time: 0:00:00  15.20 MB/s
    pysam-0.11.2.2 100% |#############################################################| Time: 0:00:00  12.58 MB/s

    (ody)$

Note the use of "conda-forge" and "bioconda".  These are "channels" that provide specialized sets of conda packages.  You can add a channel like so:

    :::bash
    $ conda config --add channels conda-forge
    $ conda config --add channels bioconda

You do not have to have a clone activated to do this; the command adds text to the `.condarc` file in your home directory.

The [bioconda channel](https://bioconda.github.io/) can be particularly useful for informatics software setup.

## Within a conda clone, you can also install packages with `pip`.
Most pip installs will be successfull:

    :::bash
    (ody) $ pip install biom-format
    Collecting biom-format
      Downloading biom-format-2.1.6.tar.gz (860kB)
        100% |████████████████████████████████| 870kB 496kB/s
    Collecting pandas>=0.19.2 (from biom-format)
      Downloading pandas-0.20.3-cp27-cp27mu-manylinux1_x86_64.whl (22.4MB)
        100% |████████████████████████████████| 22.4MB 29kB/s
    Building wheels for collected packages: biom-format
      Running setup.py bdist_wheel for biom-format ... done
    Successfully built biom-format
    Installing collected packages: pandas, biom-format
      Found existing installation: pandas 0.17.1
        Uninstalling pandas-0.17.1:
          Successfully uninstalled pandas-0.17.1
    Successfully installed biom-format-2.1.6 pandas-0.20.3
    (ody) $

## `pip` or `python setup.py` installs are needed when linking to system-tuned libraries
Many Python packages used in informatics and scientific computing in general are Python wrappers on top of high performance, system-tuned libraries like MPI or blas/lapack.  As a result, it is preferable in these cases to do a Python source install via `pip` or `python setup.py` rather than `conda install`.

For example, if you want to use MPI to parallelize some python calculations, you'll want to install mpi4py against your preferred Odyssey MPI module.

    :::bash
    $ module load python/2.7.14-fasrc01
    $ source activate ody
    (ody) $ module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
    (ody) $ pip install mpi4py
    Collecting mpi4py
    Installing collected packages: mpi4py
    Successfully installed mpi4py-2.0.0
