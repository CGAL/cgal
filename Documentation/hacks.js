function generate_autotoc() {
    var toc = $("#autotoc").append('<ul></ul>');
    if(toc) {
        toc = toc.find('ul');
        $("h1, h2, h3").each(function(i) {
            var current = $(this);
            var levelTag = current[0].tagName.charAt(1);
            current.attr('id', 'title' + i);
            toc.append("<li class='level" + levelTag + "'><a id='link" + i + "' href='#title" +
                       i + "' title='" + current.prop("tagName") + "'>" + current.text() + "</a></li>");
        });
    }
}

// hack around a bug in resize.js. The original resizeHeight is not
// considering margins or borders when calculating the div heights. Our
// ready hook is run after the original resize. We fixup the footer
// variable and then do our thing.
function fix_resize() {
    footer = $("#footer");
    resizeHeight = function() {
        var headerHeight = header.outerHeight(true);
        var footerHeight = footer.outerHeight(true);
        var windowHeight = $(window).outerHeight(true) - headerHeight - footerHeight;
        content.css({height:windowHeight + "px"});
        navtree.css({height:windowHeight + "px"});
        sidenav.css({height:windowHeight + "px",top: headerHeight+"px"});
    }
    $(window).resize(resizeHeight);
    $(window).resize();
}

// throw a stick at the modules array and hijack gotoNode 
// for our own evil purposes
$(document).ready(function() {
    if (typeof modules !== 'undefined') {
        // modules has been loaded, that means we are inside the
        // documentation of a package
        NAVTREE[0][2][1][1] = modules[0][1];
        NAVTREE[0][2][1][2] = modules[0][2];
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
    // fix_resize();
});


