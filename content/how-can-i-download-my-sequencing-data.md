Title: How can I download my sequencing data
Date: 2015-10-06 11:00
Category: FAQ
Summary:Downloading many files by right-clicking each file in a webpage can be cumbersome. You can download all your files at once with wget, scp, or cp.

Downloading many files by right-clicking each file in a webpage can be cumbersome. You can download all your files at once with the following commands<br/>

If downloading to your local computer, use wget to transfer files over HTTP:

    :::bash
    wget -r -nH --cut-dirs=1 --no-parent -e robots=off --reject="index.htm*" https://software.rc.fas.harvard.edu/ngsdata/<run_name>

If you have credentials to login to an RC host (for example, rclogin03), you can use scp to copy files to your local computer:

    :::bash
    scp -r <user_name>@rclogin03.rc.fas.harvard.edu:/n/ngsdata/<run_name> /path/to/my/directory

If transferring to a directory mounted on the cluster, the best way to transfer files is with a copy command:

    :::bash
    cp -r /n/ngsdata/<run_name> /path/to/my/directory


In all cases, replace &lt;run_name&gt; above with the name of the your sequencing run (e.g., 140523_D02345_1234_AH132DADXX), and &lt;user_name&gt; with your RC user name.