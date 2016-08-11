(function($) {

  $.fn.D3ChordChart = function( options ) {
    
    var settings = $.extend({
      status:  null
    }, options);
    
    return this.each(function() {
      
      if ($(this).data('d3chordchart')) return;
      
      var d3c  = new D3ChordChart(this,settings);
      
      $(this).data('d3chordchart',d3c);
      
      d3c.init();
      
    });
  }

  var D3ChordChart = function(ele,settings) {
    
    var obj            = this;
    
    this.ele            = ele;
    this.options        = [];

    this.dataset        = [];
    this.colors         = [];

    this.padding        = 20;
    this.w              = 920;
    this.h              = 720;
    
    this.startdate      = "1969-12-31";
    this.source         = null;

    this.init = function() {

      if (settings.hasOwnProperty('options')) {
	this.options = settings['options'];
      }

      if ($(this.ele).attr('w')) {
	this.w = $(this.ele).attr('w');
      }

      if ($(this.ele).attr('h')) {
	this.h = $(this.ele).attr('h');
      }

      if ($(this.ele).attr('startdate')) {
        this.startdate = $(this.ele).attr('startdate');
      }
      if ($(this.ele).attr('source')) {
        this.source  = $(this.ele).attr('source');
      }

      obj.id     = $(obj.ele).attr('id');

      this.update_data();
  
      window.setInterval(this.update_data,600000); 
    }

    this.update_data = function() {
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

       $(ele).html("");
       if (! ele.title) {
          ele.title = "Minilims Chord Diagram";
       } 
       $(ele).html("<h3>"+ele.title+" since " + obj.startdate + "</h3>");


       $("<div id=\"last_updated\">"+date_html+"</div>").appendTo($(ele));


       var url = "../../misc/getChordData.php?startdate="+obj.startdate;

       if (obj.source) {
          url = "/scripts/minilims/chord.php?source="+obj.source+"&startdate="+obj.startdate;
       }
       console.log("URL " + url);
      var aj = $.ajax({	
   
	url: url, //"../../misc/getChordData.php?startdate="+obj.startdate,
    	async: true,
      }).done(function(text) {
	console.log(text);
	var vals = jQuery.parseJSON(text);
	
	var types = [];
	var data  = [];

	for (var key in vals) {
	  types.push(key);
	}

	for (var i in types) {
	  var tmpvals = vals[types[i]];
	  data[i] = [];

	  for (var j in types) {
	    if (tmpvals[types[j]]) {
	      data[i][j] = parseInt(tmpvals[types[j]]);
	     // console.log("Type ",types[i]," ",types[j]," ",tmpvals[types[j]]);
	    } else {
	      data[i][j] = 0;
	    }
	  }
	}
	obj.dataset = data;
        obj.types   = types;
	obj.plot_data();
      });
	       
    }

    this.plot_data = function() {

      var chord = d3.layout.chord()
	.padding(.05)
	.sortSubgroups(d3.descending)
	.matrix(this.dataset);

      var width  = obj.w;
      var height = obj.h;

      var innerRadius = Math.min(width, height) * .31;
      var outerRadius = innerRadius * 1.1;
      
      var fill = d3.scale.ordinal()
	.domain(d3.range(4))
	.range(["#000000", "#FFDDDD", "#954444", "#F26262"]);
      
      var svg = d3.select(ele).append("svg")
	.attr("width", width)
	.attr("height", height)
	.append("g")
	.attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

      this.svg = svg;

      svg.append("g").selectAll("path")
	.data(chord.groups)
	.enter().append("path")
	.style("fill", function(d) { return fill(d.index); })
	.style("stroke", function(d) { return fill(d.index); })
	.attr("d", d3.svg.arc().innerRadius(innerRadius).outerRadius(outerRadius))
	.on("mouseover", this.fade(.05))
	.on("mouseout", this.fade(1));
 
       r1 = outerRadius*1.2;
       svg.append("g").selectAll(".arc")
            .data(chord.groups)
            .enter().append("svg:text")
            .attr("dy", ".35em")
            .attr("text-anchor", function(d) { return ((d.startAngle + d.endAngle) / 2) > Math.PI ? "end" : null; })
            .attr("transform", function(d) {
              return "rotate(" + (((d.startAngle + d.endAngle) / 2) * 180 / Math.PI - 90) + ")"
                  + "translate(" + (r1 - 15) + ")"
                  + (((d.startAngle + d.endAngle) / 2) > Math.PI ? "rotate(180)" : "");
            })
            .text(function(d) {
                return obj.types[d.index];
            });

      var ticks = svg.append("g").selectAll("g")
	.data(chord.groups)
	.enter().append("g").selectAll("g")
	.data(this.groupTicks)
	.enter().append("g")
	.attr("transform", function(d) {
	  return "rotate(" + (d.angle * 180 / Math.PI - 90) + ")"
            + "translate(" + outerRadius + ",0)";
	});
      
      ticks.append("line")
	.attr("x1", 1)
	.attr("y1", 0)
	.attr("x2", 5)
	.attr("y2", 0)
	.style("stroke", "#000");
      
      ticks.append("text")
	.attr("x", 8)
	.attr("dy", ".35em")
	.attr("transform", function(d) { return d.angle > Math.PI ? "rotate(180)translate(-16)" : null; })
	.style("text-anchor", function(d) { return d.angle > Math.PI ? "end" : null; })
	.text(function(d) { return d.label; });
      
      svg.append("g")
	.attr("class", "chord")
	.selectAll("path")
	.data(chord.chords)
	.enter().append("path")
	.attr("d", d3.svg.chord().radius(innerRadius))
	.style("fill", function(d) { return fill(d.target.index); })
	.style("opacity", 1);
	

    }

    // Returns an array of tick angles and labels, given a group.
    this.groupTicks = function(d) {
      var k = (d.endAngle - d.startAngle) / d.value;
      return d3.range(0, d.value, 100).map(function(v, i) {
	return {
	  angle: v * k + d.startAngle,
	  label: i % 5 ? null : v 
	  //label: i % 5 ? null : v / 1000 + "k"
	};
      });
    }

    // Returns an event handler for fading a given chord group.
    this.fade = function(opacity) {
      return function(g, i) {
	obj.svg.selectAll(".chord path")
          .filter(function(d) { return d.source.index != i && d.target.index != i; })
	  .transition()
          .style("opacity", opacity);
      };
    }
  }
}


)(jQuery);