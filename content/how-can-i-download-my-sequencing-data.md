Title: How can I download my sequencing data
Date: 2020-11-03 00:00
Modified: 2022-03-11
Category: FAQ
Summary:Downloading many files by right-clicking each file in a webpage can be cumbersome. You can download all your files at once with wget, scp, or cp.

Downloading many files by right-clicking each file in a webpage can be cumbersome.
The following methods are recommended (from most-to-least preferable) for downloading all your files at once.

_Excluding BCL files_

The `Data/` directory contains the raw BCL files from the sequencer.
These can take almost as much disk space as the demultiplexed FASTQ files, increasing transfer times and local storage requirements.
Each transfer mechanism downloads the entire sequencing run directory by default, but can be adapted to exclude the Data subdirectory if desired.

## Globus (with or without a FAS RC account)

[Globus](https://www.globus.org/data-transfer) provides a user-friendly web interface and the fastest and most-reliable mechanism for transfer sequencing data files within or off the FAS RC cluster.
The destination must have a Globus [endpoint](https://docs.globus.org/faq/globus-connect-endpoints/#what_is_an_endpoint) (note [destination path restrictions apply for FAS RC Holyoke and Boston endpoints](https://docs.rc.fas.harvard.edu/kb/globus-file-transfer/#Using_the_Harvard_FAS_RC_Holyoke_or_Boston_Endpoints)).
[Globus Connect Personal](https://www.globus.org/globus-connect-personal) can be used to create a Globus endpoint for download to a workstation or standalone server.

_Excluding BCL files_

Uncheck the box next to the `Data` subdirectory before clicking the "Start" button to initiate the transfer.

---

## With a FAS RC acccount

If you have a [FAS RC account](https://docs.rc.fas.harvard.edu/kb/how-do-i-get-a-research-computing-account/), several additional tools for copying data are available:

### cp

If you are transferring to another directory mounted on the FAS RC cluster, the `cp` command copy command can be used to copy the files; e.g.:

    :::bash
    cp -r /n/ngsdata/<run_name> /path/to/my/directory

### scp

If you have credentials to login to an RC host, you can use [scp or sftp](https://docs.rc.fas.harvard.edu/kb/copying-data-to-and-from-cluster-using-scp/) to copy files to your local computer; e.g.:

    :::bash
    scp -r <user_name>@login.rc.fas.harvard.edu:/n/ngsdata/<run_name> /path/to/my/local/directory

In all cases, replace &lt;run_name&gt; above with the name of the your sequencing run (e.g., 140523_D02345_1234_AH132DADXX), and &lt;user_name&gt; with your RC user name.


---

## Without a FAS RC account

If you do not have login access to the FAS RC cluster, sequencing data can be downloaded via HTTPS from https://data.rc.fas.harvard.edu/ngsdata/ via several softare tools, including:

## rclone

[rclone](https://rclone.org/) is an open-source command-line utility that can be used to transfer files over a number of storage protocols, including [HTTP directory listings](https://rclone.org/http/).
While much slower than Globus, rclone is more-performant than `wget` for transferring data from https://data.rc.fas.harvard.edu/ngsdata/ due to its ability to transfer multiple files concurrently.

    :::sh
    rclone copy --progress --http-url https://data.rc.fas.harvard.edu/ngsdata/ :http:<run_name> <destpath>

e.g., to transfer the contents oof 211007_A00794_0504_AHLMC5DSX2 locally to into a directory of the same name (creating if it doesn't already exist)

    :::sh
    rclone copy --progress --http-url https://data.rc.fas.harvard.edu/ngsdata/ :http:211007_A00794_0504_AHLMC5DSX2 211007_A00794_0504_AHLMC5DSX2

_Excluding BCL files_

Add the `--exclude='Data/**'` option:

    :::sh
    rclone copy --exclude='Data/**' --progress --http-url https://data.rc.fas.harvard.edu/ngsdata/ :http:<run_name> <destpath>

## wget

If downloading to your local computer, [wget](https://www.gnu.org/software/wget/) (installed by default on some Linux distributions, and available for [Windows](http://gnuwin32.sourceforge.net/packages/wget.htm), as well as macOS through package managers such as [conda](https://anaconda.org/conda-forge/wget) and [brew](https://formulae.brew.sh/formula/wget)) can be used to transfer files over HTTPS:

    :::bash
    wget -r -nH --cut-dirs=1 --no-parent -e robots=off  --no-check-certificate --reject="index.htm*" https://data.rc.fas.harvard.edu/ngsdata/<run_name>/

**NOTE: the trailing slash ("&lt;run_name&gt;/") is necessary to download only &lt;run_name&gt; and not the entire contents of ngsdata/**

We use InCommon certificates that may or may not be part of the trusted authorities on your local machine, so --no-check-certificate may be necessary.

_Excluding BCL files_

Add the `--exclude-directories=Data` option to the `wget` command.

