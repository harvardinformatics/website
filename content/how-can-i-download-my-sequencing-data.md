---
Title: How can I download my sequencing data
Date: 2020-11-03 00:00
Modified: 2022-03-15
Category: FAQ
Summary: List of recommended software tools for downloading sequencing data
---

[TOC]

#### Supported Protocols

Several protocols are available for downloading a sequencing run directory.
The following table contains a partial listing of software tools for each protocol.

speed | protocol | software tools | FAS RC account required
---|---|---|---
fastest | local file copy | __fpsync*__, rsync, cp | yes
| Globus | __Globus web app*__ | no
| SCP/SFTP   | __FileZilla*__, rsync, scp | yes
 slowest | HTTPS      | __rclone*__, wget | no

__*__ _recommended_

_Excluding BCL files_

The `Data/` subdirectory contains the raw BCL files from the sequencer.
These can take almost as much disk space as the demultiplexed FASTQ files, increasing transfer times and local storage requirements.
Each transfer mechanism downloads the entire sequencing run directory by default, but can be adapted to exclude the `Data/` subdirectory if desired.

# Local file copy

Users with a [FAS RC account](https://docs.rc.fas.harvard.edu/kb/how-do-i-get-a-research-computing-account/) who intend to copy the sequencing run directory (`/n/ngsdata/<run_name>`) to another file system on the FAS RC cluster may [connect to the cluster via ssh](https://docs.rc.fas.harvard.edu/kb/terminal-access/) and perform a local file copy using one of several software tools.

### rsync

_**What**_

[rsync](https://docs.rc.fas.harvard.edu/kb/rsync/) is a command-line utility for copying/synchronizing a source directory with a destination.
Unlike `cp` or `scp`, `rsync` can be efficiently resumed if transfer is interrupted.

`rsync` is installed by default on macOS.


_**Who**_

All users with a [FAS RC account](https://docs.rc.fas.harvard.edu/kb/how-do-i-get-a-research-computing-account/) who either:

* Intend to copy the sequencing run directory from `/n/ngsdata/<run_dir>` to another FAS RC Cluster file system directory, _or_
* Download to a local server / workstation, _and_ prefer a command-line utility

_**How**_

If you logeged into the FAS RC cluster via SSH, to copy files to another directory mounted on the FAS RC cluster:

    rsync -av /n/ngsdata/<run_name> /path/to/my/directory

*Note: adding a trailing slash to `<run_name>/` will cause only the contents of `<run_name>` to be copied, excluding the `<run_name>` directory itself.*

_Excluding BCL files_

Add the `--exclude=Data` option:
    
    rsync --exclude=Data -av /n/ngsdata/<run_name> /path/to/my/directory

### fpsync
[fpsync](https://www.fpart.org/fpsync/) can be used to parallelize rsync transfers (see FAS RC [Transferring Data on the Cluster](https://docs.rc.fas.harvard.edu/kb/transferring-data-on-the-cluster/) guide).
By default, fpsync uses 2 (rsync) workers to copy the directory.
The `-n <number_of_workers>` option can be used to increase the number of concurrent rsync transfers:

    fpsync -v -n 4 /n/ngsdata/<run_name> /path/to/my/directory/<run_name>

*Note: `fpsync` copies the **contents** of `<run_name>` to the destination path; specify `<run_name>` at the end of the destination path to copy the contents into a directory called `<run_name>`.*

---

# Globus

_**What**_

[Globus](https://www.globus.org/data-transfer) provides a user-friendly web interface and the fastest and most-reliable mechanism for transfer sequencing data files off the FAS RC cluster.

There is also a [Globus Command Line Interface](https://docs.globus.org/cli/) (not discussed here).

_**Who**_

All users with:

1. A Harvard Key or other supported [organizational login](https://app.globus.org/); or a Google account, ORCID iD, or [Globus ID](https://www.globusid.org/what); _and_
2. A destination Globus [endpoint](https://docs.globus.org/faq/globus-connect-endpoints/#what_is_an_endpoint); either:
    - A workstation or standalone server with [Globus Connect Personal](https://www.globus.org/globus-connect-personal) installed (available for Windows, macOS, and Linux), _or_
    - A data transfer node (maintained by your organization's system administrators) with Globus Connect Server
        - _FAS RC Cluster users_: See FAS RC [Globus File Transfer](https://docs.rc.fas.harvard.edu/kb/globus-file-transfer/) guide; note [destination path restrictions apply for FAS RC Holyoke and Boston endpoints](https://docs.rc.fas.harvard.edu/kb/globus-file-transfer/#Using_the_Harvard_FAS_RC_Holyoke_or_Boston_Endpoints).

_**How**_

1. In the demultiplex summary email, click on direct link to the run directory in the Globus web app.
   You will be prompted to authenticate using your chosen identity provider (Harvard users: select "Harvard University" in the "Look up your organization..." dropdown menu)
    - Alternatively, you may first log in to the [Globus web app](https://app.globus.org), then search for the "Harvard Bauer Core Sequencing Results" Collection, and finally double-click on the run folder that was listed in the demultiplex summary email.
2. In the other Collection box, search for your destination endpoint (which may be a Globus Connect Personal endpoint created during installation to a local workstatio/server).
   Choose or create an appropriate destination folder.
3. Select all files/folders to transfer (a checkbox in the top-left corner of the file selector can be used to select all), and select "Start" to begin the transfer.

_Excluding BCL files_

Uncheck the box next to the `Data` subdirectory before clicking the "Start" button to initiate the transfer.

---
# SSH

Users with a FAS RC account who intend to intend to transfer data off the FAS RC cluster can use one of several software tools that use the SSH protocol.

### FileZilla

_**What**_

[FileZilla](https://filezilla-project.org/) is a graphical SFTP/FTPS utility.

_**Who**_

All users with a [FAS RC account](https://docs.rc.fas.harvard.edu/kb/how-do-i-get-a-research-computing-account/) who intend to download to a local workstation.

_**How**_

1. Follow the [FAS RC SFTP file transfer using Filezilla](https://docs.rc.fas.harvard.edu/kb/sftp-file-transfer/) to install/configure FileZilla.
2. Specify `/n/ngsdata/<run_dir>` for the "Remote site"
3. Select all files/directories; drag & drop to an appropriate "Local Site" directory

_Excluding BCL files_

Unselect the `Data` subdirectory in the remote site window before dragging & dropping the selected files/folders to the local site window.


### rsync

In addition to local copies, `rsync` can also be used to push/pull data from the FAS RC cluster to a local file system or remote server that has SSH enabled.

__If you are logged into a local workstation / server, to pull data via rsync over ssh__

    rsync -av <fasrc_username>@login.rc.fas.harvard.edu:/n/ngsdata/<run_name> /path/to/my/directory

__If you are logged into the FAS RC cluster, to push data via rsync over ssh to a remote server that you have an account on__

    rsync -av /n/ngsdata/<run_name> <remote_username>@<remote_host>:/path/to/my/remote/directory

### scp

_**What**_

* The `scp` and `sftp` command-line utilities are installed by default on [Windows](https://docs.microsoft.com/en-us/windows-server/administration/openssh/openssh_overview), macOS, and most Linux distributions.

_**Who**_

* All users with a [FAS RC account](https://docs.rc.fas.harvard.edu/kb/how-do-i-get-a-research-computing-account/) who intend to download to a local workstation or another server.

_**How**_

[scp or sftp](https://docs.rc.fas.harvard.edu/kb/copying-data-to-and-from-cluster-using-scp/) can be used to a run directory to a local computer; e.g., the following `scp` command invoked from your local computer:

    scp -rp <fasrc_username>@login.rc.fas.harvard.edu:/n/ngsdata/<run_name> /path/to/my/local/directory

From an SSH session on the FAS RC cluster, to a copy files to a remote server that you have SSH access to:

    scp -rp /n/ngsdata/<run_name> <user_name>@<remote_hostname>:/path/to/my/remote/directory


In all cases, replace &lt;run_name&gt; above with the name of the your sequencing run (e.g., 140523_D02345_1234_AH132DADXX), and &lt;user_name&gt; with your RC user name.

_Excluding BCL files_

scp does not provide an option to exclude directories.


---
# HTTPS

If you do not have login access to the FAS RC cluster, sequencing data can be downloaded via HTTPS from https://data.rc.fas.harvard.edu/ngsdata/ with one of several softare tools, including:

## rclone

__**What**__

[rclone](https://rclone.org/) is an open-source command-line utility that can be used to transfer files over a number of storage protocols, including [HTTP directory listings](https://rclone.org/http/).
While much slower than Globus, rclone is more-performant than `wget` for transferring data from https://data.rc.fas.harvard.edu/ngsdata/ due to its ability to transfer multiple files concurrently.

__**Who**__

* Users who do not have login access to the FAS RC cluster, and would prefer to pull sequencing data from a web server via HTTPS; _or_
* Users who have a FAS RC account, who would like to push data from the FAS RC cluster to cloud storage (see [rclone – transfer files to/from cloud storage](https://docs.rc.fas.harvard.edu/kb/rclone/) for an example of how to configure rclone to push data to Google Drive).

__**How**__

To download https://data.rc.fas.harvard.edu/ngsdata/`<run_name>`

    rclone copy --progress --http-url https://data.rc.fas.harvard.edu/ngsdata/ :http:<run_name> <destpath>

e.g., to transfer the contents of 211007_A00794_0504_AHLMC5DSX2 locally to into a directory of the same name (creating if it doesn't already exist)

    rclone copy --progress --http-url https://data.rc.fas.harvard.edu/ngsdata/ :http:211007_A00794_0504_AHLMC5DSX2 211007_A00794_0504_AHLMC5DSX2

_Excluding BCL files_

Add the `--exclude='Data/**'` option:

    rclone copy --exclude='Data/**' --progress --http-url https://data.rc.fas.harvard.edu/ngsdata/ :http:<run_name> <destpath>

## wget

_**What**_

If downloading to your local computer, [wget](https://www.gnu.org/software/wget/) (installed by default on some Linux distributions, and available for [Windows](http://gnuwin32.sourceforge.net/packages/wget.htm), as well as macOS through package managers such as [conda](https://anaconda.org/conda-forge/wget) and [brew](https://formulae.brew.sh/formula/wget)) can be used to transfer files over HTTPS:

_**Who**_

* Users without a FAS RC account, _and_ who prefer to not use rclone (faster) or Globus (fastest).

    wget -r -nH --cut-dirs=1 --no-parent -e robots=off  --no-check-certificate --reject="index.htm*" https://data.rc.fas.harvard.edu/ngsdata/<run_name>/

**NOTE: the trailing slash (`<run_name>/`) is necessary to download only `<run_name>` and not the entire contents of `/ngsdata/**`

We use InCommon certificates that may or may not be part of the trusted authorities on your local machine, so --no-check-certificate may be necessary.

_Excluding BCL files_

Add the `--exclude-directories=Data` option to the `wget` command.