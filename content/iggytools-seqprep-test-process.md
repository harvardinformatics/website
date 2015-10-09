Title: Iggytools SeqPrep test process
Date: 2015-01-03 11:00
Category: LIMS
Summary: Procedure for testing the SeqPrep code

1. Testing seqprep needs some raw run directories which tend to be quite big so aren’t in the main github repo.  These are currently sitting on the system in <span style="line-height: 1.714285714; font-size: 1rem;"> /n/informatics/git/testdata/</span>
    
    The current runs are
    <table><tr><td>NextSeq</td><td>150527_NS500422_0126_AH2LC5AFXX</td></tr>
    <tr><td>HiSeq</td><td>150622_SN343_0517_AC6MUEACXX</td></tr>
    </table>
        
    There is also a directory for SampleSheets both correct and flawed in /n/informatics/git/testdata/SampleSheets.
    
    Make sure you’re on a machine that can see these directories.  Sandy2 is always a good choice.

1. Check out the iggy code from the git repo.   (I usually put it in /n/informatics/mclamp/test and have to do this from sa01)
    
    <pre>
    git clone git@github.com:harvardinformatics/IggyTools
    cd IggyTools
    </pre>
    
1. Setup the environment
    Ordinarily you would run the setup.sh script ( <pre>. setup.sh</pre>) if you were going to run seqprep for real.  This sets up the environment with the real environment (the real primary data and analysis directories).  We need to change these defaults before running the tests.
    
    The main config file for seqprep is usually found in
    
    <pre>~/.iggytools/iggyseq_settings.yaml</pre>
    
    There is a test version of this file in the git repo in
    
    <pre>IggyTools/tests/data/iggytools_prefs</pre>

    To make seqprep read this set an environment variable
    
    <pre>export IGGYPREFDIR=IggyTools/test/data/iggytools_prefs</pre>
    
    You’ll also need to edit the iggyseq_settings.yaml file to change the location of the following :
    
    <pre>
    PRIMARY_PARENT
    
    USERS_FILE
    
    WATCHERS_FILE
	
    LOGDIR_PARENT
    	
    PROCESSING_PARENT
    
    FINAL_PARENT
    </pre>
    
1. Copy over some test files
    
    <pre>cp -r /n/informatics/saved/ /n/informatics/mclamp/test/saved/</pre>

1. Then run
    
    <pre>python IggyTools/iggytools/bin/seqprep</pre>

    and it should give you some help text
