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

$(document).ready(function() {
    // set-up footnote generation
    $("#doc-content").append('<ol id="autoFootnotes0" class="footnotesList"></ol>');
    $("body").footnotes();
    generate_autotoc();
    // fix_resize();
});
