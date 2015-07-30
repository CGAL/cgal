function generate_autotoc() {
    var toc = $("#autotoc").append('<ul></ul>');
    if(toc.length > 0) { // an autotoc has been requested
        toc = toc.find('ul');
        var indices = new Array();
        indices[0] = 0;
        indices[1] = 0;
        indices[2] = 0;

        $("h1, h2, h3").each(function(i) {
            var current = $(this);
            var levelTag = current[0].tagName.charAt(1);
            var cur_id = current.attr("id");

            indices[levelTag-1]+=1;  
            var prefix=indices[0];
            if (levelTag >1) prefix+="."+indices[1];
            if (levelTag >2) prefix+="."+indices[2];
            current.html(prefix + "   " + current.html());
            for(var l = levelTag; l < 3; ++l){
                indices[l] = 0;
            }

            if(cur_id == undefined) {
                current.attr('id', 'title' + i);
                current.addClass('anchor');
                toc.append("<li class='level" + levelTag + "'><a id='link" + i + "' href='#title" +
                           i + "' title='" + current.prop("tagName") + "'>" + current.text() + "</a></li>");
            } else {
                toc.append("<li class='level" + levelTag + "'><a id='" + cur_id + "' href='#title" +
                           i + "' title='" + current.prop("tagName") + "'>" + current.text() + "</a></li>");
            }

        });
    }
}

// throw a stick at the modules array and hijack gotoNode 
// for our own evil purposes
$(document).ready(function() {
    if (typeof modules !== 'undefined') {
        // modules has been loaded, that means we are inside the
        // documentation of a package
        NAVTREE[0][2][1][1] = modules[0][1];
        NAVTREE[0][2][1][2] = modules[0][2];
        // override gotoNode from navtree.js
        gotoNode = function (o,subIndex,root,hash,relpath) {
            var nti = navTreeSubIndices[subIndex][root+hash];
            if(nti && (nti[0] === 1 && nti[0])) {
                nti.splice(1, 1);
            }
            o.breadcrumbs = $.extend(true, [], nti ? nti : navTreeSubIndices[subIndex][root]);
            if (!o.breadcrumbs && root!=NAVTREE[0][1]) { // fallback: show index
                navTo(o,NAVTREE[0][1],"",relpath);
                $('.item').removeClass('selected');
                $('.item').removeAttr('id');
            }
            if (o.breadcrumbs) {
                o.breadcrumbs.unshift(0); // add 0 for root node
                showNode(o, o.node, 0, hash);
            }
        }
    }
    // set-up footnote generation
    $("#doc-content").append('<ol id="autoFootnotes0" class="footnotesList"></ol>');
    $("body").footnotes();
    generate_autotoc();
});


/* 
 * A jQuery plugin by Brian Holt that will search the selected blocks for 
 * specially-defined footnote elements.	 If found, these elements will be 
 * moved to a footnotes section and links to and from the footnotes will 
 * be created.
 *
 * See http://www.planetholt.com/articles/jQuery-Footnotes	
 * for full documentation.
 *
 * By default, footnotes will be found in SPANs with the footnote class,
 * and in BLOCKQUOTEs with a TITLE attribute.
 *
 * Thanks to CSSNewbies.com for the general idea, which I have enhanced 
 * and implemented with as a jQuery plugin.
 * 
 * Copyright 2008-2009 Brian Holt.	
 * Licensed under the LGPL license. See 
 * http://www.gnu.org/licenses/lgpl-3.0-standalone.html
 * 
 * Version 1.2.2
 */
(function(c){c.fn.footnotes=function(d){var e=c.extend({},c.fn.footnotes.defaults,d);return this.each(function(f){b("INFO: Building footnotes for "+(f+1)+"...",e.debugMode);c(e.footnotes,this).addClass(e.autoFootnoteClass);var h=(""===e.contentBlock)?c(this):c(e.contentBlock,this),g=e.orderedList?"<ol/>":"<ul/>";c("."+e.autoFootnoteClass).each(function(l){var t=-1,n=f+"-"+l,q=c(this),j,r,s,u,p,m,o,k;if(e.singleFootnoteDestination){j=c("#"+e.destination);if(0===j.length){b("INFO: No #autoFootnotes found; adding our own",e.debugMode);j=c(g).attr("id",e.destination).addClass("footnotesList").appendTo(h)}}else{j=c("#"+e.destination+f);if(0===j.length){b("INFO: No #autoFootnotes"+f+" found; adding our own for "+(f+1),e.debugMode);j=c(g).attr("id",e.destination+f).addClass("footnotesList").appendTo(h)}}q.removeClass(e.autoFootnoteClass);r=e.fnExtractFootnote(this);t=-1;n=f+"-"+l;j.find("li > .footnoteContent").each(function(i){var v=c(this);if(v.html()===r){t=i;s=c(v.parents("li").get(0));return false}});if(-1===t){u=c("<a/>").attr("href","#cite-text-"+n).attr("name","cite-ref-"+n).attr("id","cite-ref-"+n).attr("dir","ltr").attr("title",r).text("["+(j.find("li").length+1)+"]").addClass("footnoteLink");if(q.is(e.prependTags)){c("<sup/>").prependTo(this).append(u)}else{c("<sup/>").appendTo(this).append(u)}p=c("<li/>").attr("id","cite-text-"+n);m=c("<span/>").addClass("footnoteBackReferenceGroup").appendTo(p);c("<span/>").addClass("footnoteContent").html(r).appendTo(p);u=c("<a/>").text("^").attr("href","#cite-ref-"+n).addClass("footnoteBackref").prependTo(m);j.append(p)}else{n=f+"-"+t;m=c(c("li > .footnoteBackReferenceGroup",j).get(t));o=m.find(".footnoteBackref");k=o.length;if(0===o.length){b("ERROR: $backRefs.length == 0, which should have prevented this code path",e.debugMode)}else{if(1===o.length){c("<sup/>").text("^ ").addClass("footnoteBackref").prependTo(m);o.html("<sup>a</sup>");++k}u=c("<a/>").attr("href","#"+s.attr("id")).attr("name","cite-ref-"+n+"-"+o.length).attr("id","cite-ref-"+n+"-"+o.length).attr("title",r).text("["+(t+1)+"]").addClass("footnoteLink");if(q.is(e.prependTags)){c("<sup/>").prependTo(this).append(u)}else{c("<sup/>").appendTo(this).append(u)}u=c("<a/>").attr("href","#cite-ref-"+n+"-"+o.length).addClass("footnoteBackref");if(k>=26){b("WARN: multiple letter functionality is probably broken when more than 26 footnotes exist",e.debugMode)}u.prepend(String.fromCharCode((k)+96));c("<sup/>").appendTo(m).append(u)}}});b("INFO: Done building footnotes for "+(f+1),e.debugMode)})};c.fn.footnotes.version=function(){return"1.2.2"};c.fn.footnotes.defaults={footnotes:"blockquote[title],span.footnote,blockquote[cite]",prependTags:"blockquote",singleFootnoteDestination:false,destination:"autoFootnotes",contentBlock:".content",autoFootnoteClass:"autoFootnote",fnExtractFootnote:a,orderedList:true,debugMode:true};function b(e,d){if(d){if(window.console&&window.console.log){window.console.log(e)}}}function a(i){var j=c(i),e,f,h,g,d;if(j.is("span.footnote")){e=c(i).html();f=/^(?:(?:&nbsp;)|\s)*\(([\S\s]+)\)(?:(?:&nbsp;)|\s)*$/;h=e.match(f);if(h&&2===h.length){e=h[1]}j.empty()}else{if(j.is("blockquote[title]")){g=j.attr("cite");e=j.attr("title");if(""!==g){d=c("<a/>").attr("href",g);if(0===c(e).length){e=d.text(e)}else{e=d.text(g).wrap("<span/>").parent().append(": "+e);j.attr("title","")}}}else{if(j.is("blockquote[cite]")){g=j.attr("cite");e=c("<a/>").attr("href",g).text(g)}}}return e}})(jQuery);
