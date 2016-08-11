(function($){
	$(document).ready(function(){
		$("#toc-title").click(function(){
			if ($("#toc-toggle").html() == "hide") {
				$("#toc-toggle").html("show");
				$("#toc").hide();
			}
			else {
				$("#toc-toggle").html("hide");
				$("#toc").show();
			}
		});
	});
})(jQuery)