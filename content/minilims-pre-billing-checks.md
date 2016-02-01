Title: Minilims - Pre-Billing Checks
Date: 2016-01-01
Category: Tutorials
Tags: MiniLIMS
Author: Michele Clamp
Summary: Billing checks to do before invoicing.

Before invoices can be generated there are a set of billing checks that can be run that check everything is valid and can proceed to invoicing.  These things are : 

* Each group has the right status (HARVARD,OUTSIDE_ACADEMIC,COMMERCIAL) 
* Each Submission has valid billing information whether it is expense code(s) or Purchase Orders.  

If the Submission is split across expense codes then it checks everything adds up to 100%.

    :::bash
    php plugins/Billing/InvoiceGeneratorTests.php -m mm -y yyyy

If errors are found it will report things like

### Missing Group Type

    No group type for group [Wirth] Allowed values are COMMERCIAL,HARVARD,OUTSIDE_ACADEMIC

To correct this - go to the relevant group page (click groups on sidebar and then the group name),  press the edit button and select the right status.   If it's not obvious then have a look at the affiliation of some of the group members to work it out.

### Wrongly formatted Expense Codes - fixable

    Setting submission Expense_Code_1 from [NSF EAGER : 385-34240-6600-144199-363225-0001-10005] to [NSF EAGER : 385-34240-8250-144199-363225-0001-10005]

These corrections (changing 6600 to 8250 or adding in dashes to the expense code) won't be saved by default.   Use the -s option with the php command to save the changes.

### Wrongly formatted Expense Codes - not fixable

Other expense code booboos include

    [0] => Wrong expense code format for [Reagent_Request][REA00803] - [370.31660.825.133190.340437.0001.45042]
    [1] => Expense code doesn't match pattern [370.31660.825.133190.340437.0001.45042] for sub [REA00803]
    [2] => Wrong expense code format for [Reagent_Request][REA00828] - [370316506600000770600200000045023]
    [3] => Expense code doesn't match pattern [370316506600000770600200000045023] for sub [REA00828]

For these you have to go into the Submission (or in this case the Reagent_Request) and reformat them manually.

### Missing Purchase Orders

    [0] => No purchase order [0006707982] exists for [Reagent_Request][REA00796]

For these you have to create a new Purchase_Order entry (<span style="color: #ff0000;">page</span>) if one doesn't already exist.
