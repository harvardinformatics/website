/*
 * Minilims StatusBox
 *
 * When applied to a div it will query a minilims server with a type and property name and return
 * a json structure with all values of the property and their count 
 * 
 */

(function($) {


  /* The plugin is structured with two main parts
   *  
   * 1. The html element    (var element)  - contains the gui components (div)
   * 2. The StatusBox object (var obj)     - contains the functions and the data
   *
   * Usage:
   * 
   * To create :
   * $("#statusbox").MinilimsStatusBox();
   *
   */

  $.fn.MinilimsStatusBox= function(options) {

    return this.each(function() {

      if ($(this).data('statusbox')) return;

      var statusbox = new MinilimsStatusBox(this);

      $(this).data('statusbox',statusbox);

      statusbox.init();

    });
  }

  var MinilimsStatusBox = function(element) {
    var obj     = this;
    var ele     = element;

    this.init = function() {
       obj.count     = 0;
       this.type     = $(element).attr('type');
       this.property = $(element).attr('property');
       this.source   = $(element).attr('source');
       this.title    = $(element).attr('title');
       this.color    = $(element).attr('color');
       this.baseurl  = "/scripts/minilims/status.php";

       this.updateData();

       window.setInterval(this.updateData,60000); 
    }

    this.updateData = function() {

       var host = null;

       var url =  obj.baseurl + "?";
       url += "source="   +obj.source;
       url += "&type="    +obj.type;
       url += "&property="+obj.property;
       url += "&startdate=previousmonth";

       console.log(url);

       var data;
  
       var d = new Date();
       var curr_date = d.getDate();
       var curr_month = d.getMonth();
       var curr_year = d.getFullYear();
       var curr_hour = d.getHours();
       var curr_minute = d.getMinutes();
       if (curr_minute < 10) {
	   curr_minute = "0" + curr_minute;
       }
       var date_html = "<h5>Last updated: " + curr_month + "/" + curr_date + "/" + curr_year + " " + curr_hour + ":" + curr_minute + "</h5>";

       $("#last_updated").html(date_html);

       $(ele).html("<h3>"+obj.title+"</h3>");
       $.ajax({
          url:  url,
          data:  data,
          success:  function(data) {
               console.log(data);
               var plotdata = [];

               var str = "<h1>"+obj.title+"</h1>";
               var i = 0;

               $.each( data, function( key, val ) {
                   var tmp = {};

                   if (val > 0) {
                   tmp['key'] = key;
                   tmp['val'] = val;
                   tmp['i']   = i; 
                   i++;
                     console.log("counting " + key);
                   } else {
                     console.log("Not counting " + key);
                   }
                   plotdata.push(tmp);
               });

               console.log(plotdata);

               var margin = {top: 10, right: 10, bottom: 80, left: 20};
	       var width  = 150 - margin.left - margin.right;
	       var height = 200 - margin.top  - margin.bottom;

	       var x = d3.scale.ordinal().rangeRoundBands([0, width], .05);
	       var y = d3.scale.linear().range([height, 0]);

	       var xAxis = d3.svg.axis()
		   .scale(x)
		   .orient("bottom");
		   
	       var yAxis = d3.svg.axis()
		   .scale(y)
		   .orient("left")
		   .ticks(5);

	       var svg = d3.select(ele).append("svg")
		   .attr("width", width + margin.left + margin.right)
		   .attr("height", height + margin.top + margin.bottom)
		   .append("g")
		   .attr("transform",
			 "translate(" + margin.left + "," + margin.top + ")");

	       x.domain(plotdata.map(function(d) { return d.key; }));
	       y.domain([0, d3.max(plotdata, function(d) { return d.val; })]);

	       svg.append("g")
		   .attr("class", "x axis")
		   .attr("transform", "translate(0," + height + ")")
		   .call(xAxis)
		   .selectAll("text")
		   .style("text-anchor", "end")
		   .attr("dx", "-.8em")
		   .attr("dy", "-.55em")
		   .attr("transform", "rotate(-90)" );

	       svg.append("g")
		   .attr("class", "y axis")
		   .call(yAxis)
		   .append("text");

	       svg.selectAll("bar")
		   .data(plotdata)
		   .enter().append("rect")
		   .style("fill", obj.color)
		   .attr("x", function(d) { return x(d.key); })
		   .attr("width", x.rangeBand())
		   .attr("y", function(d) { return y(d.val); })
		   .attr("height", function(d) { return height - y(d.val); });
	       
	       },
          async: true,
          dataType: "json"
      });
      obj.count++; 
  }

 }
})(jQuery);
