Title: Adding a new user email to an existing alert
Date: 2015-03-12 00:00
Category: Tutorials
Tags: Helium
Author: Michele Clamp
Summary: How to add a new user to a Helium MiniLIMS alert.

Minilims can be set up to send email when new pages are added or pages change. To see all the existing alerts click the 'Alerts' link on the left hand sidebar.


<figure>
	    <a class="img" href="/images/helium-email1.png">
		    <img class="img-responsive" src="/images/helium-email1.png"></img>
	    </a>
<figcaption></figcaption>
</figure>   

This should show you a list of alerts in a table.

<figure>
	    <a class="img" href="/images/helium-email2.png">
		    <img class="img-responsive" src="/images/helium-email2.png"></img>
	    </a>
<figcaption></figcaption>
</figure>   

Select the Alert you want to change. For instance if you want to add a user to be notified whenever a new Dewar is requested then select the alert with Page_Type = Helium_Dewar_Request, Trigger_Property = Status and Trigger_Value = REQUESTED.

In this case the Alert is ALE00001

Once on the Alert page click the Edit button and you'll get a form where you can edit the fields.

Add a new Email address (separated by a comma from previous addresses) into the'Action_Value' text box. Click the Save button to save and you are all set.

<figure>
	    <a class="img" href="/images/helium-email3.png">
		    <img class="img-responsive" src="/images/helium-email3.png"></img>
	    </a>
<figcaption></figcaption>
</figure>   
