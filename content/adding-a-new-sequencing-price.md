Title: Adding a new Sequencing Price
Date: 2016-01-01
Category: Tutorials
Tags: Next-Gen Sequencing, Bauer Core, MiniLIMS
Author: Michele Clamp
Summary: How to add a new Sequencing Price to the Bauer Core MiniLIMS system


## Create a new Sequencing Price record

[Log in](https://bauer-minilims.rc.fas.harvard.edu/minilims "Log in") as an admin into minilims.

Go to the 'Add Page->HarvardSCLims->Sequencing Price' menu entry to bring up a price form. 

<figure>
	<a class="img" href="/images/sequencing-price-1.png">
    		<img class="img-responsive" src="/images/sequencing-price-1.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Fill out the form with the relevant fields and prices.  If the status is set to inactive then it won't be used for any billing or cost estimations. 

<figure>
	<a class="img" href="/images/sequencing-price-2.png">
    		<img class="img-responsive" src="/images/sequencing-price-2.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

The Display_Name field has to be in a certain format for everything to be billed correctly.  We have to match the run type selected in the submission form to this field in each price.

## Construct the Display_Name Field

For new prices the format of the Display_Name should be

    <run_type> <readcount> x <readlen> bp

* `<run_type>` is the sample run type in the [submission form](https://bauer-minilims.rc.fas.harvard.edu/minilims//plugins/Core/submit_batch.php?parenttype=Submission&childtype=Sample&name=new "submission form") (see screenshot below)
 
<figure>
	<a class="img" href="/images/sequencing-price-3.png">
    		<img class="img-responsive" src="/images/sequencing-price-3.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

* readcount is 1 for unpaired or 2 for paired - read length should be a number e.g.

    Megaseq (v2) 2 x 100 bp

Megaseq (v2) is the run type selected in the submission form 2 means a paired end run 100 is the read length. Please make sure there is only one space between the elements. The different run_types are found in the submission form (and can be configured from the Configure Page->Sample configuration screen)

## Existing Display_Name Fields are slightly different

For historical reasons some of the sequencing prices are treated differently as the format is a little different. These are listed below:

### NextSeq Runs

    Submission run type      Price Display_Name

    NextSeq_High             "NextSeq <1|2> x <readlen> bp - High"
    NextSeq_Mid              "NextSeq <1|2> x <readlen> bp - Mid"

### Rapid Runs

    Submission run type     Price Display_Name

    Rapid_single_lane       “Rapid Run <1|2> x <readlen> bp - Single Lane"
    Rapid_full_flow_cell    “Rapid Run <1|2> x <readlen> bp - Full Flowcell"

### Standard Runs

    Submission run type    Price Display_Name
    Standard_v4            “Standard Run (v4) <1|2> x <readlen> bp"
    Standard_v3            “Standard Run <1|2> x <readlen> bp"