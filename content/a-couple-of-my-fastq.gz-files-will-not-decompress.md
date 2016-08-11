Title: A couple of my fastq.gz files will not decompress, or give errors when I process them. What is wrong?
Date: 2015-10-06 11:00
Category: FAQ
Summary:Occasionally, large files become corrupted during download. You can determine whether this has occurred by comparing the <a href="http://en.wikipedia.org/wiki/Checksum">checksum</a> of a file before and after download. We have placed checksums for your fastq.gz files in your run directory in a file called md5sum.txt. Compare the values in this file to new checksums calculated on your downloaded files.

Occasionally, large files become corrupted during download. You can determine whether this has occurred by comparing the [checksum](http://en.wikipedia.org/wiki/Checksum) of a file before and after download. We have placed checksums for your fastq.gz files in your run directory in a file called md5sum.txt. Compare the values in this file to new checksums calculated on your downloaded files.

To calculate a checksum for a file called myfile.fastq.gz, use:

    :::bash
    md5sum myfile.fastq.gz
    
If the checksums for a file before and after download differ, then the file is corrupt or incomplete and should be downloaded again.

