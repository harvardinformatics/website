Title: Basic Unix Workshop
Date: 2016-01-13 13:51
Category: Tutorials
Author: mclamp@harvard.edu
Summary: This tutorial is intended for people completely new to unix and command line driven analysis.  It introduces how to navigate the unix filesystem, make and delete files and directories and be able to list the contents of files and directories.


## Prerequisites

*   A Harvard FAS RC cluster account 
    [https://account.rc.fas.harvard.edu/request/](account-request>)

*   A terminal program so you can log into the cluster 
    [https://rc.fas.harvard.edu/docs/access-and-login.html](access-and-login>)

*   The ability to log into the cluster 
    [https://rc.fas.harvard.edu/docs/access-and-login.html](access-and-login>)



These are also covered in the first 14 slides of the RC Intro to Unix slides here

[https://software.rc.fas.harvard.edu/training/intro_unix/latest/#(1)](https://software.rc.fas.harvard.edu/training/intro_unix/latest/#(1))

This is a very good intro put together by John Brunelle and is worth a read. We’ll be covering most of the things but with more of an emphasis on file content manipulation.

## Summary of Commands Covered

* pwd 
* cd 
* mkdir
* cp
* ls
* man
* rmdir
* rm 
* mv
* cat
* head
* tail
* less
* scp

## 1. Logging in and the command prompt

When you first successfully log into the cluster you will see a window containing something similar to

    :::bash
    [mclamp@rclogin04 ~]$

The mclamp will be replaced by your username and the rclogin04might be different as there are a number of login nodes available and you'll land on one of them at random. This string is called the 'command prompt' and when you see it you are able to type in commands.

The `~` is shorthand for your home directory. If you move out of your home directory then this will be replaced by the name of the directory you're working in.

The region after the $ is where you type your commands. Any output from the command is displayed on your terminal window and when the command completes the command prompt is displayed again ready for the next one. 

For long-running commands you have the ability to run these in the
background and have the command prompt returned to you straight away. For this tutorial all commands should only take a second or so to process so we won't need this ability.

<div class="exercises">
Exercises
<ul>
<li> Log into the cluster now. How does your command prompt differ from the above
example? Does this make sense?</li>
<li> Press the return key at your command line - what happens?</li>
<li> Type a simple command whoami followed by the return key at your command line.</li>
What happens? Is this true?
</ul>
</div>

## 2. Directories

### Directory Structure

Like windows and macs unix filesystems are divided into directories in a hierarchical fashion. Directories are separated by a forward slash (/) and all files hang off a single top-level directory called root (written as /).

### Home Directory

When you first log in you land in your home directory. This has a full path like `/n/home01/msmith`. When you want to use your home directory in a command you can either type in the full path or you can refer to it via the shortcut `~`.


### Working Directory

As you move around the filesystem your working directory will change. Your command prompt will give you the last 'chunk' of your working directory but the pwd (print working directory) command will give you the full path. For example :

    :::bash
    [mclamp@rclogin04 ~]$ pwd
    /n/home_rc/mclamp
    [mclamp@rclogin04 ~]$


<div class="exercises">
Exercises
<ul><li>What is the full path of your home directory?</li></ul>
</div>

## 3. Useful Filesystems

Even though all directories hang off the root / directory there are many different types of filesystems on the cluster. Some are backed up and have hourly snapshots whereas others are not backed up but are very fast and suitable for intermediate processing results.

### Home directories 

backed up, have snapshots in the .snapshot directory, relatively slow and
shouldn't be used for large scale cluster usage. 40Gb quota.

### Scratch directories. 

These are very fast but not backed up filesystems. There is no charge for usage but people need to email rchelp@fas.harvard.edu to request access. There is no limit to usage but files unused for 3 months will be purged from the system. The main scratch filesystem is `/n/regal`

### Lab shares. 

These are backed up and can be mounted on your laptop/desktop for easy access to files. These are usually not suitable for large scale analysis and we strongly recommend people use the scratch systems for processing and then copy over output files to their lab shares for safekeeping.

## 4. Command line editing

Unix work involves a fair bit of typing which, although powerful, is prone to typing mistakes. The bash shell has a number of built-in utilities that help reduce errors and speed up your work.

### Tab Completion

If you start typing a command or a filename and then press TAB twice the shell will print out to the terminal a list of possible completions. This works with both commands and with files. For example if you type the letters pw and then TAB TAB you will get a list of possible command completions.

    :::bash
    [mclamp@rclogin04 ~]$ pw
    pwck pwconv pwd pwdx pwhich pwmconfig pwunconv
    [mclamp@rclogin04 ~]$ pw


Here we can see there are a number of possible commands that start with the letters pw. If there is only one possible command the shell fills in the whole command for you and doesn't present a list of options. e.g. if we type whoa TAB TAB we get

    :::bash
    [mclamp@rclogin04 ~]$ whoami


This also works very well for files and directories. If you start typing a filename and then press the tab key twice you'll get a list of possible completions. As with the command completion if there is only one possibility it will fill the whole filename in for you.


For instance if you type `ls .bash` TAB TAB in your home directory (ls is the command to list files)
you get:

    :::bash
    [mclamp@rclogin04 ~]$ ls .bash
    .bash_history .bash_logout .bash_profile .bashrc
    [mclamp@rclogin04 ~]$ ls .bash


<div class="exercises">
Exercises
<ul>
<li>Using the ls command try tab completion in the root directory and your home directory</li>
<li>Do you need to even start typing a file to use tab complete (try ls SPACE TAB TAB to
see)</li>
<li>How many commands starting ls are there on the system?</li>
</ul>
</div>

### Line Navigation

The unix command line is text based and you can't use the mouse to click and position the cursor. Instead you have to use the keyboard. At first this may well seem cumbersome but there are only a few really useful keystrokes to learn and once learned become extremely fast and useful. The most commonly used strokes are:

<table style="margin: 1em">
<tr><td style="width: 200px">Left and right arrow keys</td><td>Move the cursor left and right</tr>
<tr><td><code>Ctrl-d</code></td><td>Delete character under cursor</td></tr>
<tr><td><code>Ctrl-a</code></td><td>Move to the beginning of the line</td></tr>
<tr><td><code>Ctrl-e</code></td><td>Move to the end of the line</td></tr>
<tr><td>Up arrow</td><td>Previous command (you can go through the whole command history this way)</td></tr>
<tr><td>Down arrow</td><td>Next command</td></tr>
<tr><td><code>Ctrl-k</code></td><td>Delete from cursor to the end of the line</td></tr>
</table>


Still useful but less commonly used keystrokes

<table style="margin: 1em">
<tr><td><code>Ctrl-w</code></td><td>Delete to beginning of word</td></tr>
<tr><td><code>Ctrl-r</code></td><td>Search for most recent command containing</td></tr>
<tr><td><code>Ctrl-xx</code></td><td>Move between current cursor position and beginning of line</td></tr>
</table>

You can combine these together. For instance a useful one is `Ctrl-a Ctrl-k` which goes to the beginning of the line and deletes it.

<div class="exercises">
Exercises
<ul>
<li>Cut and paste this command to your command line ls -ltra /tmp/ |sort -n ></li>
/tmp/tmp.dat Use this for each of the following exercises and try to complete in as
few keystrokes as you can.
<li>Change -ltra to -s</li>
<li>Change /tmp/ to /usr</li>
<li>Change sort -n to wc -l</li>
<li>Replace everything after ls with ~</li>
<li>Repeat the command starting ls -ltra /usr/</li>
</ul>
</div>

### Command Control Shortcuts

Up to now we've only used commands that return a small amount of output and run very fast. Quite often this won't be the case and a command will start to print pages of output to the screen or take a long time to run (or both). In these cases there are some commands that can halt the output to the screen so you can read it and also to kill the command completely.

<table style="margin: 1em">
<tr><td><code>Ctrl-s</code></td><td>Stop the output to the screen</td></tr>
<tr><td><code>Ctrl-q</code></td><td>Continue the output to the screen</td></tr>
<tr><td><code>Ctrl-c</code></td><td>Kill the currently running command</td></tr>
<tr><td><code>Ctrl-z</code></td><td>Suspend the currently running command (follow by bg to push into the background)</td></tr>
<tr><td><code>Ctrl-l</code></td><td>Clears the screen (Useful if your terminal fills up and becomes confusing</td></tr>
</table>

<div class="exercises">
Exercises
<ul>
<li>Enter the command ls -l /usr/bin and then clear the screen</li>
<li>Run the top command which shows you the top running processes. Enter the keystroke</li>
<li> kill it and return to the command line</li>
</ul>
</div>

### Command History

Another feature of the bash shell to save typing is a set of keystrokes/commands to locate and rerun
previous run commands. These include


* `history` Prints out all stored previous commands
* `!!` Run last command
* `!mycmd` Run the most recent command that starts with mycmd (e.g. !cd)
* `!mycmd:p` Print out the command that !mycmd would run (also adds it as the latest command in the command history)

<div class="exercises">
Exercises
<li>Enter the command to list all stored previous commands</li>
<li>Rerun the previous command with as few keystrokes as possible</li>
<li>Rerun the most recent command that starts with ls</li>
<li>Rerun the most recent command that starts with ls but change the second argument to</li>
your home directory
</ul>
</div>

## 5. Getting Help

You can find out how to use every unix command by typing man (man stands for manual). These man pages give a brief description of the function of each command and then usually a list of options and how to use them. This is the first part of the ls man page for listing files:

    :::man
    LS(1) User Commands LS(1)
    NAME
     ls - list directory contents
    SYNOPSIS
     ls [OPTION]... [FILE]...
    DESCRIPTION
     List information about the FILEs (the current directory by default). Sort entries
     alphabetically if none of -cftuvSUX nor --sort.
     Mandatory arguments to long options are mandatory for short options too.
     -a, --all
     do not ignore entries starting with .
     -A, --almost-all
     do not list implied . and ..
     --author
     with -l, print the author of each file
     -b, --escape
     print octal escapes for nongraphic characters


Most man pages don't fit on one screen and you can navigate through them using the following
commands

* `spacebar` Show next page
* `Ctrl-b` Show previous page
* Down/up arrow Scroll down/up one line
* `/mystr` Search for mystr (press return after the string to search)
* `q` Quit the man page display and go back to the command prompt

<div class="exercises">
Exercises
<ul>
<li>Enter the command to show the man page for the lscommand</li>
<li>Use the scrolling commands to find the option to list one file per line</li>
<li>Use the search function to find the option to list without sortinge</li>
</ul>
</div>


## 6. Directory Navigation

Now we come to some more commands. First we need to be able to move about the filesystem. The main command for this is `cd` (change directory). It is used as follows:


    :::bash
    cd /n/regal/
    cd ~
    cd /n/regal/informatics_public/ref

The tab completion is very useful here. You can just start a directory name (say `/n/ref`) and press the tab key twice and it'll either complete it for you (if it's unique) or provide you with a list of options.

The `cd` command will always get you where you need to go but there are two other useful commands that are useful when you just need to go to a directory quickly and then back again. There are `pushd` and `popd`. They are used as follows:

    :::bash
    [mclamp@sandy2 ~]$ pushd /tmp/
    /tmp ~
    
    [mclamp@sandy2 tmp]$ ls
    76738.unp hsperfdata_mclamp pip-build-root
    hsperfdata_afreedman krb5cc_57100_O6SyG1 sess_39g827l9ecj3evijrih06t2cu6
    [mclamp@sandy2 tmp]$ popd
    ~
    [mclamp@sandy2 ~]$

This lets you go back to your original directory without typing the whole name in. For your home directory it's not a hardship as you can use ~but if you're deep in a directory tree it can save a lot of typing.

You don't have to use full paths for directory navigation. If you want to move relative to your current directory the `..` notation will take you one level up. For example :

    :::bash
    cd /n/regal/informatics_public/ref/ensembl
    cd ../ucsc

is the same as

    :::bash
    cd /n/regal/informatics_public/ref/ensembl
    cd /n/regal/informatics_public/ref/ucsc

Similarly the single dot `.` notation refers to the current directory. This is useful when copying files

<div class="exercises">
Exercises
<ul>
<li>Use the cd command to move to /n/regal/informatics_public/ref/ensembl/directory.</li>
<li>Use the pushd/popd commands to move to your home directory, use the ls command and then move back to the original directory</li>
<li>Check that the .. notation really does move you one level up</li>
</ul>
</div>


## 7. Listing files

We've used this command in previous sections but now we'll go into a little more detail. This is a command you'll use a lot for listing files in a directory. It has a lot of options to change the display and show different kinds of information about files.

To just show a list of filenames in a directory do the following

    :::bash
    ls <dirname>

For example

    :::bash
    ls /n/regal
    [mclamp@sandy2 ~]$ ls /n/regal/
    Giribet_lab     bertoldi_lab    chien_lab       desai_lab       edwards_lab
    Meissner_lab    betley_lab      clardy_lab      desmarais_lab   eggan_lab
    TSC             bomblies_lab    combes_lab      dsouza_lab      eisenstein_lab
    adams_lab       branton_lab     conroy_lab      dulac_lab       ellison_lab
    airoldi_lab     brenner_lab     coull_lab       duraisingh_lab  engert_users
    ascherio_lab    capellini_lab   cox_lab	        dutton_lab      ervin_lab
    baccarelli_lab  chetty_lab      davis_lab       economics       evans_lab
    berger_lab      chevrier_lab    demler_lab      eddy_lab        extavour_lab

This is a bare bones listing and doesn't distinguish between files and directories. We can add options to the `ls` command to show different things

    :::bash
    ls -F /n/regal
    Giribet_lab/    bertoldi_lab/   chien_lab/      desai_lab/      edwards_lab/
    Meissner_lab/   betley_lab/     clardy_lab/     desmarais_lab/  eggan_lab/
    TSC/            bomblies_lab/   combes_lab/     dsouza_lab/     eisenstein_lab/
    adams_lab/      branton_lab/    conroy_lab/     dulac_lab/      ellison_lab/
    airoldi_lab/    brenner_lab/    coull_lab/      duraisingh_lab/ engert_users/
    ascherio_lab/   capellini_lab/  cox_lab/        dutton_lab/     ervin_lab/
    baccarelli_lab/ chetty_lab/     davis_lab/      economics/      evans_lab/
    berger_lab/     chevrier_lab/   demler_lab/     eddy_lab/       extavour_lab/

Here directories are shown with a trailing slash to distinguish them from files. The most important option to know about for the ls command is -l which shows the long listing. This is an example:

    :::bash
    [mclamp@sandy2 ~]$ ls -l /n/informatics/mclamp/iggy/
    total 164
    drwxrwxr-x 8 mclamp rc_admin 327 Jun 13 12:32 IggyTools
    -rw-rw-r-- 1 mclamp rc_admin 146 Aug 10 16:29 README.md
    drwxrwxr-x 5 mclamp rc_admin 67 May 12 13:28 iggy_ve2.7.3
    -rwxrwxr-x 1 mclamp rc_admin 314 Aug 17 12:08 setup.example
    -rwxrwxr-x 1 mclamp rc_admin 314 Aug 17 12:08 setup.sh

Here we have multiple columns for each file/directory. Below are descriptions of what they contain

    -rw-rw-r-- 1 mclamp rc_admin 146 Aug 10 16:29 README.md
    <----1---> 2 <------3------> <4> <----5-----> <----6--->


####  1. Permissions  (`drwx-rw-r--`)
This defines the permissions for the file or directory. There are 10 characters and they represent
different types of permission

* `drwx-rw-r--` The first character is usually whether it is a directory (d) or a plain file (-)
* `drwxrw-r--` The next 3 characters define the permissions for the owner r = can read, w = can write, x = can execute
* `drwxrw-r--` The next 3 characters similarly define permissions for the group
* `drwxrw-r--` The next 3 characters similarly define permissions for everyone


#### 2. Inodes  (`1`)
The next number is the number of inodes used for the file (google this if you're interested)

#### 3. Ownership (`mclamp rc_admin`)
The next two columns are the owner's username and the group

#### 4. Size (`146`)
The next column is the size of the file in bytes

#### 5. Last time accessed (`Aug 10 16.29`)
The next 3 columns are the time the file was last accessed

#### 6. Filename (`README.md`)
The final column is the filename itself.

### Listing a subset of files

With no extra options the ls command returns a list of everything. You can return subsets by using partial filenames and wildcards. For example:

* `ls -l *.fastq` lists all files with the .fastq extension (* is the wildcard)
* `ls -l B*.txt` lists all files beginning with B and having the .txt extension

### Listing all files

Unix has the concept of 'hidden' files that aren't listed by default. These are all files/directories starting with `..` For instance in your home directory you likely have a .bashrc file which stores configuration for your shell. You can list everything in a directory using the -a option

* `ls -la` List all files including hidden ones


### Sorting a file list
By default `ls -l` returns a list in alphabetical order. Often you want to know the oldest/newest file or the largest/smallest file. Common sort options are:

* `ls -lt` List all files sorted by modification time
* `ls -lX` List all files sorted by extension
* `ls -lS` Sort by file size
* `ls -lrS` Reverse sort by file size (-r reverses all sorts


### Human readable format for file size
File sizes can get very large and listing the size in bytes can get confusing. Using the -h (for human) makes the listing more readable and the size will be formatted into 100 (bytes) 100k (kbytes) 100M (megabytpes) 100G (gigabytpes) 100T (terabytes).

* `ls -lh` List files with sizes in human readable format

<div class="exercises">
Exercises
<ul>
<li>How many directories are in the listing for /n/regal/informatics/workshops/Basic_UNIX_2015-11-04?</li>
<li>What is the largest file in /n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data/seq? </li>
<li>Using the man page find the option to ls that recursively lists files in a directory.</li>
<li>Which is the biggest file in the whole directory tree under /n/regal/informatics/workshops/Basic_UNIX_2015-11-04?</li>
</ul>
</div>

## 8. Creating files

### 8.1 Using an editor

For newcomers to the cluster we recommend the editor nano to create and edit files. Just type `nano` on the command line to start up the program.

Handily many of the useful commands within nano are listed at the bottom of the screen (The ^ character means the control key).

### 8.2 Using the cat command

Using an editor is usually the way to go when you want to create a text file but you can use the cat command with the redirect character > to quickly add text to a file.

    :::bash
    cat > myfile.dat
    Add some text in here
    And some more on a new line if you like.
    To exit press ctrl-d
    ^D


Once you've created it you can again use the cat command to check the content is what you wanted it to be

    :::bash
    cat myfile.dat


### 8.3 Using redirect from a command into a file

Very often we are not entering in text from the keyboard but want to capture output from a command. We will go into more detail about redirects at the end of the workshop but the simplest form is

    :::bash
    mycommand > myfile.dat

For example

    :::bash
    ls -l ~ > myhomedir.txt


This will put the output from the ls command into a file myhomedir.txt

<div class="exercises">
Exercises
<ul>
<li>Using the cat command and then redirecting to a file create a file myscript.sh in your
home directory containing the contents of /n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data/myscript.sh</li>
<li>Using the nano editor edit the file you just made in your home directory. This is an
example of a shell script. This contains a series of commands (the top line #!/bin/bash
says use the bash shell) that are run when the file is invoked. Edit the content of the file
to add 2 extra commands 

    <ul><li>1) Change to the /n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data directory</li>
        <li> 2) Do a long listing of that directory</li>
    </ul>

<li>Run the script you just made (bash ./myscript.sh) and redirect the output to
myscript.out in your home directory. Look at the output - is it what you expect? </li>
</ul>
</div>

## 9. Looking at File Contents

We used one of the ways to look at a file's contents in the previous section - the cat command. This is fine for short files that fit into the window but we're often dealing with large files, sometimes many gigabytes. For these occasions one of the following options is better.

### 9.1 The head command - beginnings of files

The head command will by default just show you the top 10 lines of any file

    :::bsh
    head myscript.sh

If we want to see more lines we can give it the number of lines to show with the `-n` option

    :::bash
    head -n 100 myfile.dat

This will show the first 100 lines of myfile.dat

At first glance this may not seem particularly useful but in practice is a very handy way of just peeking into a file to see what the contents are. It is also useful for checking that the output from a command is producing what you expect

### 9.2 The tail command - ends of files

Similarly the tail command shows by default the last 10 lines of a file. And again you can use the `-n` option to display more lines

    :::bash
    tail -20 myfile.dat


One very handy option to tail is the -f option. This will show the last 10 lines of a file but if the file is appended to will also show those lines as well. This is extremely useful to show the growth of an output file over time.

    :::bash
    tail -f mygrowingfile.dat

You can type ctrl-c to exit the command

### 9.3 The less command - scrolling and searching large files

The cat, head and tail commmands are useful but sometimes you just want to look at a file, scroll up and down and be able to search. In these situations less is the command you need and it has many options and is very powerful. We are only going to cover a few of the most useful ones.

The basic syntax for less is

    :::bash
    less mybigfile.dat


This will show you the first page of mybigfile.dat and you can navigate using the following commands


<br>
<ul>
<li>SPACE Scroll down half a page</li>
<li>up/down arrows Scroll up/down 1 line</li>
<li>q Quit</li>
<li>h Shows the help page</li>
<li>ctrl-D Scroll down half a page</li>
<li>ctrl-U Scroll up half a page</li>
<li>G Go to end of file</li>
<li>gg Go to beginning of file</li>
<li>ctrl-g Show status line at bottom of screen</li>
<li>10g Go to line 10</li>
</ul>


In addition to navigating files we can search for patterns


<br>
<ul>
<li>/mystring Goes to the line with the next occurence of mystring</li>
<li>!mystring Goes to the line *not* containing mystring</li>
<li>n Finds the next occurence of the search string</li>
<li>N Finds the previous occurence of the search string</li>
<li>&mystring Shows only those lines containing mystring (& followed by return to get out)</li>
<li>&!mystring Shows only those lines *not* containing mystring</li>
</ul>

<div class="exercises">
Exercises
<li>Use the tail command to see if the file
/n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data/mylongscript.out is
changing and roughly how often</li>
<li>Look at the source for this file /n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data/mylongscript.sh . Use the man command to find out what this output means</li>
<li>Practice the less navigation commands on the file in
/n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data/seqruns.dat</li>
<li>Use the navigation commands to find what is on the first and last lines</li>
<li>Find the first occurence of the string 1501 (Jan 2015). Roughly how far through the file
is this</li>
<li>Use the search commands to only show lines with the number of Read 1 Cycles. What is
the largest read length?</li>
<li>How many Kb were sequenced in the 10th run in the file?</li>
</ul>
</div>


## 10. Making Directories - `mkdir`

The `mkdir` command will create a new directory. For example

    :::bash
    mkdir ~/mynewdir

will create a new directory in your home directory.

If you want to create a directory multiple levels down use the -p option to fill in the middle directories.

    :::bash
    mkdir -p ~/dir1/dir2/dir3


<div class="exercises">
Exercises
<ul>
<li>Make yourself a user directory (use your username as the directory name) in</li>
<li>/regal/informatics/workshops/Basic_Unix_2015-11-04/Users</li>
<li>Use a single command to create under the directory you just created called
test1/test2/test3</li>
</li>
</div>


## 11. Copying Files and Directories - `cp`

Continuing the tradition of terse commands copying files and directories is done using the `cp` command. For files the basic syntax is :

    :::bash
    cp myfile.dat newfile.dat

We can also copy files into directories using

    :::bash
    cp myfile.dat mydir


Note that the directory has to exist for this command to work.

Multiple files can be copied using wildcards.

    :::bash
    cp *.dat mydir


This will copy all files ending in .dat into a new directory.

For directory copying we generally want to copy the contents as well. In this case we use the -r option

    :::bash
    cp -r mydir mynewdir


In this case the destination directory doesn't need to pre exist.

By default the `cp` command makes 'fresh' copies of all the files and doesn't preserve timestamps. If you want to do this use the `-p` option.

    :::bash
    cp -rp mydir mynewdir


Finally sometimes we want to reference a file in a directory but don't want to copy the whole thing. This is useful to save space for very large files and for this we use the `ln` command to link to the original.

    :::basj
    ln -s myrealfile.dat mylinkedfile.dat


If you then do an `ls` to list the linked file you'll see something like

    :::bash
    [mclamp@bioinf02 ~]$ cd work
    [mclamp@bioinf02 work]$ ln -s ../data/Homo_sapiens.GRCh38.78.gtf human.gtf
    [mclamp@bioinf02 work]$ ls -l *.gtf
    lrwx------ 1 mclamp rc_admin 34 Nov 2 16:08 human.gtf -> ../data/Homo_sapiens.GRCh38.78.gtf*

<div class="exercises">

Exercises
<ul>

<li>Make a Data directory under your username directory you created in the last section.</li>
<li>Find the largest file in the /n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data directory</li>
<li>Make a symlink to this file in your newly created Data directory</li>
<li>Use head or less to make sure you have the right content in the file and your symlink worked</li>
</ul>
</div>


## 12. Deleting Files and Directories - `rm`

To remove a file use the `rm` command

    :::bash
    rm myfile.dat


Note if you're on a scratch filesystem you won't be able to get this back. If you're on a backed up filesystem you can retrieve things from the .snapshot directories.

To remove a directory use the `rm -rf` command which recursively removes the directory and all it's contents.

    :::bash
    rm -rf myunwanteddir

Note - be very careful with this and double check before pressing return. Be especially careful with wildcards.

<div class="exercises">
Exercises
<ul>
<li>Copy the /n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data/seq
directory and its contents to your user space</li>
<li>Change into that directory and remove all the files starting with pog</li>
<li>Remove your seq directory and the contents with a single command</li>
</ul>
</div>

## 13. Renaming Files and Directories - `mv` and `rename`

You can rename files and directories using `mv`

    :::bash
    mv myfile.dat mynewfile.dat
    mv mydir mynewdir

You can also use the relative notation to move things up a directory or directories

    :::bash
    mv myfile.dat ../../

This will move myfile.dat up two directories

If you want to rename multiple files then you can't just use wild cards. You either have to use a loop or the rename command. An example of the rename command is :

    :::bash
    rename .htm .html *.htm

This will rename all files matching *.htm from an .htm extension to an .html extension. Using a loop is a little more complicated but the following will do the same thing

    :::bash
    for i in *.htm ; do # Loop over all files ending .htm and reference them using $i
     j=${i%.htm} # Strip off the .htm extension and reference that using $j
     mv $i $j.html # Use the mv command to rename the original file to .html
    done

<div class="exercises">
Exercises
<ul>
<li>Copy the /n/regal/informatics/workshops/Basic_UNIX_2015-11-04/Data/seq
directory and its contents to your user space</li>
<li>Change to that directory and rename all the .fa files to .fasta files</li>
<li>Now change them all back again</li>
</ul>
</div>


## 14. Transferring Files To/From the Cluster

There are graphical tools to transfer files to and from the cluster but it is very handy to know the command line versions. For this we use the `scp` (secure copy) command.

The basic syntax is (and you need a terminal open on your local machine):

    :::bash
    scp username@server:remotefilename localfilename

The remote filename can be the full path or just the filename if it's in your home directory. The local filename can just be . for your current directory.

For example for a single file we can use:

    :::bash
    scp mclamp@login.rc.fas.harvard.edu:/n/mypath/seq/pog.fa .

For a directory and its contents we can use the -r option :

    :::bash
    scp -r mclamp@login.rc.fas.harvard.edu:/n/mypath/seq/ .

Common mistakes include

* Forgetting the colon : separating the servername from the path
* Mistyping the path
* Forgetting the -r when copying directories
* Running scp from the remote machine and not your local machine - check your command prompt people! Everyone does it at least once.

<div class="exercises">
Exercises
<ul>
<li>Transfer the seq directory used above to your local machine. Try and get it right first time</li>
</ul>
</div>
