Title: PICRUSt on Odyssey
Author: Aaron Kitzmiller
Date: 2017-07-10 00:00
Category: Software
Tags: PICRUSt, Python, Odyssey, Anaconda
Summary: Setup PICRUSt on Odyssey

From the [Picrust website](http://picrust.github.io/picrust/):

    PICRUSt (pronounced “pie crust”) is a bioinformatics software package 
    designed to predict metagenome functional content from marker gene (e.g., 16S rRNA) 
    surveys and full genomes.

## Installation
Because PICRUSt is a Python package, it can be [setup on Odyssey using a local Anaconda clone](Python on Odyssey).  It is a bit more complex than a simple pip install, though.

### Setup a conda clone 
Follow the instructions [for setting up an Anaconda clone on Odyssey](anaconda-on-odyssey), if you don't already have one. 
{% include 'anaconda-setup.md' %}

### Setup h5py
As noted in the [instructions for setting up h5py on Odyssey](h5py-on-odyssey), h5py can be pip installed into your conda clone, but you should choose a module for HDF5 based on your likely compiler and MPI library stack.
{% include 'h5py-setup.md' %}

### Install biom-format
This can be installed with a simple pip install with your clone activated

    :::bash
    (ody) $ pip install biom-format==2.1.5

### Install PICRUSt
This must be installed from an downloaded archive.  The archive need not be unpacked.  PyCogent will also be installed.

    :::bash
    (ody) $ wget https://github.com/picrust/picrust/releases/download/1.1.1/picrust-1.1.1.tar.gz
    --2017-07-10 13:57:57--  https://github.com/picrust/picrust/releases/download/1.1.1/picrust-1.1.1.tar.gz
    Resolving github.com... 192.30.253.112, 192.30.253.113
    Connecting to github.com|192.30.253.112|:443... connected.
    ---
    Saving to: “picrust-1.1.1.tar.gz”
    2017-07-10 13:57:58 (17.0 MB/s) - “picrust-1.1.1.tar.gz” saved [10826010/10826010]

    (ody) $ pip install picrust-1.1.1.tar.gz 
    Processing ./picrust-1.1.1.tar.gz
    Collecting cogent==1.5.3 (from PICRUSt==1.1.1)
      Downloading cogent-1.5.3.tgz (3.4MB)
        100% |████████████████████████████████| 3.4MB 182kB/s 
    Building wheels for collected packages: PICRUSt, cogent
      Running setup.py bdist_wheel for PICRUSt ... done
      Running setup.py bdist_wheel for cogent ... done
    Successfully built PICRUSt cogent
    Installing collected packages: cogent, PICRUSt
      Found existing installation: cogent 1.9
        Uninstalling cogent-1.9:
          Successfully uninstalled cogent-1.9
    Successfully installed PICRUSt-1.1.1 cogent-1.5.3
    
    (ody) $ format_tree_and_trait_table.py --help
    Usage: format_tree_and_trait_table.py [options] {-t/--input_tree INPUT_TREE -i/--input_trait_table INPUT_TRAIT_TABLE}

    [] indicates optional input (order unimportant)
    {} indicates required input (order unimportant)

