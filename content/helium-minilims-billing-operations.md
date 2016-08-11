Title: Helium MiniLIMS - Billing Operations
Date: 2015-03-12 00:00
Category: Tutorials
Tags: Helium
Author: Michele Clamp
Summary: Billing through the Helium MiniLIMS


The Helium LIMS tracks dewars (Helium_Dewar, Helium_Dewar_Request) that have been requested and delivered.  When they are delivered they are weighed and  the reading saved Helium_Dewar_Reading).  At the end of the month these delivered dewars are charged by the volume of Helium delivered.  This is calculated by converting from weight to volume of Helium liquid. 

A number of emergency dewars may also have been taken from the LISE building and need to be charged.  Shanna Losee will email me with details of those weights although it is very simple to enter those in by hand (<span style="color: #ff0000;">page</span>). 

Finally we need to track how much helium has been recovered.   Each group has a recovery node (or nodes) (Helium_Node) and the meter reading is uploaded at the beginning of each month (Helium_Node_Reading) (<span style="color: #ff0000;">page</span>). 


The Invoices are generated using the Invoice Generator link in the left hand sidebar.  This process is identical across all lims (<span style="color: #ff0000;">page</span>).

## Monthly Process

1\.  Upload the Helium_Node_Readings. An email is sent to me and informatics@fas.harvard.edu at the beginning of every month with the readings for every recovery node.   The email comes from RENO@harvard.edu and has a subject like

<div class="codehilite"><pre>HARVARD UNIVERSITY ACS, Trend Interval Report, <span class="il">HE</span>.MONTH.METER.REPOR<wbr>T, 2/1/2015 09:15 AM</pre></div>

a) Download the report from email and scp it to helium.rc.fas.harvard.edu 

b) ssh into helium.rc.fas.harvard.edu 

c) run the reading uploader script in test mode

    :::bash
    php plugins/Helium/scripts/load_MeterReport.php -f <reportfile>

This will print out various things to the screen if it is running fine. d) run the reading uploader script and save to the database

    :::bash
    php plugins/Helium/scripts/load_MeterReport.php -f <reportfile> -s

I generally then have a quick look at the [Helium_Node_Reading summary](http://helium.rc.fas.harvard.edu/minilims-dev//plugins/Core/Templates/Default_Summary.php?type=Helium_Node_Reading) page to check things look ok. 

e) run the post-processing script to update the monthly recovered amounts for each reading. The meter report only records the total amount that has travelled through the meter and we need to record how much has been recovered per month.   As a wrinkle sometimes the meters reset themselves to 0 which results in a negative amount recovered.   A second script calculates the monthly totals taking this into account (actually fixing the negative numbers isn't in there yet)

    :::bash
    php plugins/Helium/scripts/calc_recovered_amounts.php

If all looks good then add the -s option to save.

    :::bash
    php plugins/Helium/scripts/calc_recovered_amounts.php -s

A healthy looking set of readings looks like this 

<figure>
	<a class="img" href="/images/helium1.png">
    		<img class="img-responsive" src="/images/helium1.png"></img>
	</a>
    <figcaption>Figure 1</figcaption>
</figure>


f) Enter Emergency Dewars For a Linde dewar you need to enter a dewar request and set the status to DELIVERED. For each Harvard emergency dewar you need to enter a Helium_Dewar_Request and also a Helium_Dewar_Reading to record the weight.  We need both in order to bill correctly. The emergency dewar info is usually emailed in the format


**January Helium Emergency Dewars:**

1.  1/18/15, Lukin Lab, Alp Sipahigil, dewar HUHe 100L-14, weight: 217.8 lbs
2.  1/31/15, Kim lab, Gil-Ho Lee, dewar HUHe100L-8, weight: 223.4 lbs



## Entering a request for an emergency dewar

* Click the 'Request Dewar' link in the left hand sidebar. 
* Remember the request number (something like REQ001234) which is in the top field. 
* Enter the Actual_Delivery_Date 
* Enter the Group name (Lukin).   This will autocomplete and also query for any expense codes. 
* Enter the Requested_By person.  This will also autocomplete 
* Enter the Requested_Delivery_Date (same as the actual delivery date).   The popup will be greyed out but you can enter the date by hand by typing yyyy-mm-dd in the field. 
* Select the expense code to use (if none is there this can be fixed when generating invoices) 
* Select the correct size. 
* Enter something like 'Emergency dewar entered by <yourname>' in the comments field. 
* If this is a Linde (external) dewar set the status to DELIVERED.  Dewars won't be billed unless they have this status. 
* Press save. You now need to enter the weight reading. 
* Type in the the dewar name which is something like HuHe100L14 (case is important and Shanna generally has a hypen before the last number which is not needed). This is what a successful request looks like after saving.

<figure>
	<a class="img" href="/images/helium2.png">
    		<img class="img-responsive" src="/images/helium2.png"></img>
	</a>
    <figcaption></figcaption>
</figure>
          

* Now we enter the weight reading.  Press 'Scan Dewar' on the left hand sidebar. 
* Type in the the dewar name which is something like HuHe100L14 (case is important and Shanna generally has a hyphen before the last number which is not needed). 
* Press the 'Record Weight' button. 


<figure>
	<a class="img" href="/images/helium3.png">
    		<img class="img-responsive" src="/images/helium3.png"></img>
	</a>
    <figcaption></figcaption>
</figure>
                

* Enter the dewar request number,  the date the weight was read (the delivery date), set the status to DEWAR_DELIVERED_WEIGHT and enter the actual weight.  (Do not press the update button - this tries to automatically read the weight from the scale and we're not attached to it). 

A successful form looks like this 

<figure>
	<a class="img" href="/images/helium4.png">
    		<img class="img-responsive" src="/images/helium4.png"></img>
	</a>
    <figcaption></figcaption>
</figure>

Press the 'Save' button to save.
