Title: Genome Repository Status
Date: 2016-01-01
Author: Michele Clamp
Category: Software
Tags: Genome Databases
Summary: Current status of public genome files

<script type="text/javascript">
<![CDATA[
$.get("/scripts/get_status.php", function(data) {
	  console.log( data );
	  $("#gr_status").html(data);
	});
]]>
</script>

<div id="gr_status"></div>
