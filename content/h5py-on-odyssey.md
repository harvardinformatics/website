Title: h5py on Odyssey
Author: Aaron Kitzmiller
Date: 2018-02-07 00:00
Category: Software
Tags: h5py, HDF5, Python, Odyssey, Anaconda
Summary: Setup h5py on Odyssey

## h5py from conda can be useful
[h5py](http://docs.h5py.org/en/latest/) is a package that provides a Pythonic interface for HDF5 files through the HDF5 C library.  The Anaconda distribution that we commonly use on the Odyssey cluster can provide h5py with a simple `conda install h5py`.  This will install h5py along with an hdf5 *.so library as an all-in-one solution.  

As of today, this library is based on the 1.10 version of HDF5.

    {% include 'anaconda-setup.md' %}
    
    (ody)$ conda install h5py
	Fetching package metadata ...............
	Solving package specifications: .

	Package plan for installation in environment /n/home_rc/akitzmiller/.conda/envs/ody:

	The following NEW packages will be INSTALLED:

	    h5py: 2.7.1-py27_2 conda-forge
	    hdf5: 1.10.1-0     conda-forge

	Proceed ([y]/n)? y

	hdf5-1.10.1-0. 100% |####################################| Time: 0:00:00   9.95 MB/s
	h5py-2.7.1-py2 100% |####################################| Time: 0:00:00  18.20 MB/s
 
If you are using an older Odyssey python module (e.g. python/2.7.6-fasrc01), you may get an older, 1.8 series HDF5 library


## For builds needing Odyssey HDF5 libraries, use a pip install with `--no-binary`
In some cases, you may need to build h5py against Odyssey HDF5 libraries.  There may be a different version requirement or you may need MPI-enabled HDF5.  In this case, you should avoid conda-based binary installs and use pip for a source code installation.  

pip normally installs h5py via a binary (wheel), similar to conda.  You can force a source code installation using the `--no-binary` option.  If you first load an Odyssey HDF5 module, a source code installation will find the headers and libaries that are exposed by the module

    :::bash
    $ module load gcc/7.1.0-fasrc01 mvapich2/2.3b-fasrc01 hdf5/1.10.1-fasrc01
    $ module load python/2.7.13-fasrc01
    $ source activate ody

    (ody) $ pip install --no-binary=h5py h5py

	Installing collected packages: h5py
  		Running setup.py install for h5py ... done
	Successfully installed h5py-2.7.1

	(ody) $ ldd ~/.conda/envs/ody/lib/python2.7/site-packages/h5py/h5.so

	libhdf5.so.101 => /n/sw/fasrcsw/apps/MPI/gcc/7.1.0-fasrc01/mvapich2/2.3b-fasrc01/hdf5/1.10.1-fasrc01/lib/libhdf5.so.101 (0x00002acb0ad5c000)
	libhdf5_hl.so.100 => /n/sw/fasrcsw/apps/MPI/gcc/7.1.0-fasrc01/mvapich2/2.3b-fasrc01/hdf5/1.10.1-fasrc01/lib/libhdf5_hl.so.100 (0x00002acb0b33e000)


`ldd` shows that the h5.so library that h5py uses references the Odyssey libhdf5.so.101 library.

