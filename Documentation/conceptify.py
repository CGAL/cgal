#!/usr/bin/python3

import os
import re
import glob
from pyquery import PyQuery as pq
from lxml import etree

def conceptify_nested_classes(d):
    # change nested classes to nested concepts
    nested_classes=d('a[name=nested-classes]').parent()
    nested_classes.text('Concepts')
    nested_classes=nested_classes.parent().parent().siblings()

    # we cannot use a proper selector here because PyQuery doesn't deal
    # with empty pseudo-classes
    nested_classes=nested_classes.filter(lambda i: pq(this).attr('class') == 'memitem:')
    nested_classes.children(".memItemLeft").text("concept")


def conceptify(d):
    # fix the title
    title = d(".title")
    title.html(re.sub("Class Reference", "Concept Reference", title.html()))
    # remove the include
    include_statement = d(".textblock").next()
    # should check that this is really the div we think it is
    include_statement.empty()

    conceptify_nested_classes(d)
    
def write_out_html(d, fn):
    f = open(fn, 'w')
    # this is the normal doxygen doctype, which is thrown away by pyquery
    f.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n')
    print(d, file=f)
    f.close()

# remove duplicate files
def clean_doc():
    duplicate_files=glob.glob('./output/CGAL.CGAL.*/html/jquery.js')
    duplicate_files.extend(glob.glob('./output/CGAL.CGAL.*/html/dynsections.js'))
    for fn in duplicate_files:
        os.remove(fn)

# from http://stackoverflow.com/a/1597755/105672
def re_replace_in_file(pat, s_after, fname):
    # first, see if the pattern is even in the file.
    with open(fname) as f:
        if not any(re.search(pat, line) for line in f):
            return # pattern does not occur in file so we are done.

    # pattern is in the file, so perform replace operation.
    with open(fname) as f:
        out_fname = fname + ".tmp"
        out = open(out_fname, "w")
        for line in f:
            out.write(re.sub(pat, s_after, line))
        out.close()
        os.rename(out_fname, fname)

class_files=glob.glob('./output/CGAL.CGAL.*/html/class*.html')
for fn in class_files:
    d = pq(filename=fn, parser='html')
    ident = d('#CGALConcept')
    if ident.size() == 1:
        conceptify(d);
        # in case of a second pass don't process this again
        d.remove("#CGALConcept")
    # there is a doxygen bug that prevents the correct linkage of the CGAL breadcrumb 
    ident = d('#nav-path .navelem').eq(0).children().eq(0)
    if ident and ident.attr('href') == 'namespaceCGAL.html':
        ident.attr('href', '../../CGAL.CGAL/html/namespaceCGAL.html')
    write_out_html(d, fn)

# in a group we only need to change the nested-classes
group_files=glob.glob('./output/CGAL.CGAL.*/html/group*Concepts*.html')
for fn in group_files:
    d = pq(filename=fn, parser='html')
    conceptify_nested_classes(d)
    write_out_html(d, fn)    

# fix up Files
files_files=glob.glob('./output/CGAL.CGAL.*/html/files.html')
for fn in files_files:
    d = pq(filename=fn, parser='html')
    table = d("table.directory")
    row_id=table("td.entry").filter(lambda i: pq(this).text() == 'Concepts').parent().attr('id')
    if row_id != None:
        # figure out the rowid and then drop everything from the table that matches
        table("tr").filter(lambda i: re.match(row_id + '*', pq(this).attr('id'))).remove()
        write_out_html(d, fn)

filesjs_files=glob.glob('./output/CGAL.CGAL.*/html/files.js')
for fn in filesjs_files:
    re_replace_in_file('^.*\[ "Concepts",.*$', '', fn)
    
# external is placed by doxygen to mark a class from a tagfile, this
# is more confusing then helpful in our case
re_replace_in_file('\[external\]', '', './output/CGAL.CGAL/html/annotated.html')

# fix class/concept mismatch in generated pages
relationship_pages=glob.glob('./output/CGAL.CGAL.*/html/hasModels.html')
relationship_pages.extend(glob.glob('./output/CGAL.CGAL.*/html/generalizes.html'))
relationship_pages.extend(glob.glob('./output/CGAL.CGAL.*/html/refines.html'))
for fn in relationship_pages:
    d = pq(filename=fn, parser='html')
    dts=d(".textblock .reflist dt")
    # no contents() on pyquery, do it the hard way
    dts.each(lambda i: pq(this).html(re.sub("Class ", "Concept ", pq(this).html())))
    write_out_html(d, fn)

#throw out nav-sync
all_pages=glob.glob('./output/CGAL.CGAL*/html/*.html')
for fn in all_pages:
    d = pq(filename=fn, parser='html')
    d('#nav-sync').remove()
    write_out_html(d, fn)

clean_doc()
