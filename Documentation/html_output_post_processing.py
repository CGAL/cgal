#!/usr/bin/env python2
# Copyright (c) 2012 GeometryFactory (France). All rights reserved.
# All rights reserved.
#
# This file is part of CGAL (www.cgal.org).
# You can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# Licensees holding a valid commercial license may use this file in
# accordance with the commercial license agreement provided with the software.
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# $URL:
# $Id:
#
#
# Author(s)     : Philipp Moeller

import argparse
import glob
import codecs
from lxml import etree
import os
from os import path
from pyquery import PyQuery as pq
import re
import shutil
from sys import argv
from sys import stderr

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
    title.html(re.sub("((Class)|(Struct))( Template)? Reference", "Concept Reference", title.html()))
    # remove the include
    include_statement = d(".contents").children().eq(0)
    # should check that this is really the div we think it is
    include_statement.empty()

    conceptify_nested_classes(d)

# just an alias
def conceptify_ns(d):
    conceptify_nested_classes(d)

def write_out_html(d, fn):
    f = codecs.open(fn, 'w', encoding='utf-8')
    # this is the normal doxygen doctype, which is thrown away by pyquery
    f.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n')
    f.write(d.html())
    f.write('\n')
    f.close()

def package_glob(target):
    return filter(lambda x: not './CGAL/' in x, glob.glob(target))

# remove duplicate files
def clean_doc():
    duplicate_files=package_glob('./*/html/jquery.js')
    duplicate_files.extend(package_glob('./*/html/dynsections.js'))
    duplicate_files.extend(package_glob('./*/html/resize.js'))
    duplicate_files.extend(package_glob('./*/html/stylesheet.css'))
    # kill _all_, including the one in CGAL tabs.css files
    duplicate_files.extend(glob.glob('./*/html/tabs.css'))
    # left-over by doxygen?
    duplicate_files.extend(package_glob('./*/bib2xhtml.pl'))
    duplicate_files.extend(package_glob('./*/cgal_manual.bib'))
    duplicate_files.extend(package_glob('./*/citelist.doc'))
    duplicate_files.extend(package_glob('./*/doxygen.bst'))
    duplicate_files.extend(package_glob('./*/geom.bib'))
    duplicate_files.extend(package_glob('./*/ftv2cl.png'))
    duplicate_files.extend(package_glob('./*/ftv2ns.png'))
    
    for fn in duplicate_files:
        os.remove(fn)

# from http://stackoverflow.com/a/1597755/105672
def re_replace_in_file(pat, s_after, fname):
    # first, see if the pattern is even in the file.
    with codecs.open(fname, encoding='utf-8') as f:
        if not any(re.search(pat, line) for line in f):
            return # pattern does not occur in file so we are done.

    # pattern is in the file, so perform replace operation.
    with codecs.open(fname, encoding='utf-8') as f:
        out_fname = fname + ".tmp"
        out = codecs.open(out_fname, "w", encoding='utf-8')
        for line in f:
            out.write(re.sub(pat, s_after, line))
        out.close()
        os.rename(out_fname, fname)

def is_concept_file(filename):
  if not path.exists(filename):
    return False;
  d = pq(filename=filename, parser='html')
  ident = d('#CGALConcept')
  return ident.size() == 1

def rearrange_img(i, dir_name):
    img = pq(this)
    if img.attr("src") == "ftv2cl.png":
        parser=pq(this).parent()
        for link_class in ['a.el', 'a.elRef']:
            links=parser(link_class)
            if links.size()>0 and is_concept_file(path.join(dir_name, pq(links[0]).attr("href"))):
                img.attr("src","ftv2cpt.png")
    srcpath=img.attr("src")
    img.attr("src", "../../CGAL/html/" + srcpath.split('/')[-1])

class figure_anchor_info:
  def __init__(self):
    self.next_index=1
    self.anchor_map={}


def collect_figure_anchors(i,infos):
  anchor_name=pq(this).attr('id')
  if re.match("fig__.+",anchor_name) != None:
    infos.anchor_map[anchor_name]=infos.next_index
    infos.next_index+=1

def update_figure_ref(i,infos):
  link = pq(this)
  link_name=link.text()
  if re.match("fig__.+",link_name) != None:
    link.text( "Figure "+str(infos.anchor_map[link_name]) )

def automagically_number_figure(filename):
  infos=figure_anchor_info()
  d = pq(filename=filename, parser='html')
  d('a.anchor').each( lambda i: collect_figure_anchors(i,infos) )
  write_out_html(d, filename)
  #reference manual pages might also contain references to figures
  all_pages=glob.glob(path.join(path.dirname( path.abspath(filename) ),'*.html'))
  for fname in all_pages:
    d = pq(filename=fname, parser='html')
    d('a.el').each( lambda i: update_figure_ref(i,infos) )    
    write_out_html(d, fname)

def main():
    parser = argparse.ArgumentParser(
        description='''This script makes adjustments to the doxygen output. 
It replaces some text in specifically marked classes with the appropriate text for a concept, 
removes some unneeded files, and performs minor repair on some glitches.''')
    parser.add_argument('--output', metavar='/path/to/doxygen/output', default="output")
    parser.add_argument('--resources', metavar='/path/to/cgal/Documentation/resources')
    
    args = parser.parse_args()
    resources_absdir=args.resources
    os.chdir(args.output)

    # number figure
    main_pages=package_glob('./*/html/index.html')
    for fn in main_pages:
      automagically_number_figure(fn)

    #replace icons with CGAL colored ones
    shutil.copy(path.join(resources_absdir,"ftv2cl.png"),path.join("CGAL/html/", "ftv2cl.png"))
    shutil.copy(path.join(resources_absdir,"ftv2ns.png"),path.join("CGAL/html/", "ftv2ns.png"))
    shutil.copy(path.join(resources_absdir,"ftv2cpt.png"),path.join("CGAL/html/", "ftv2cpt.png"))

    annotated_files=package_glob('./*/html/annotated.html')
    for fn in annotated_files:
      dir_name=path.dirname(fn)
      d = pq(filename=fn, parser='html')
      tr_tags = d('table.directory tr img')
      tr_tags.each(lambda i: rearrange_img(i, dir_name))
      write_out_html(d,fn)
    
    class_files=package_glob('./*/html/class*.html')
    class_files.extend(package_glob('./*/html/struct*.html'))
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
            ident.attr('href', '../../CGAL/html/namespaceCGAL.html')
        write_out_html(d, fn)

    namespace_files=package_glob('./*/html/namespace*.html')
    for fn in namespace_files:
        d = pq(filename=fn, parser='html')
        ident = d('#CGALConceptNS')
        if ident.size() == 1:
            conceptify_ns(d);
            d.remove("#CGALConceptNS")
        write_out_html(d, fn)

    # in a group we only need to change the nested-classes
    group_files=package_glob('./*/html/group*Concepts*.html')
    for fn in group_files:
        d = pq(filename=fn, parser='html')
        conceptify_nested_classes(d)
        write_out_html(d, fn)

    # fix up Files
    files_files=package_glob('./*/html/files.html')
    for fn in files_files:
        d = pq(filename=fn, parser='html')
        table = d("table.directory")
        row_id=table("td.entry").filter(lambda i: pq(this).text() == 'Concepts').parent().attr('id')
        if row_id != None:
            # figure out the rowid and then drop everything from the table that matches
            table("tr").filter(lambda i: re.match(row_id + '*', pq(this).attr('id'))).remove()
            write_out_html(d, fn)

    filesjs_files=package_glob('./*/html/files.js')
    for fn in filesjs_files:
        re_replace_in_file('^.*\[ "Concepts",.*$', '', fn)

    #Rewrite the path of some images
    re_replace_in_file("'src','ftv2",
                       "'src','../../CGAL/html/ftv2",
                       'CGAL/html/dynsections.js')

    # external is placed by doxygen to mark a class from a tagfile, this
    # is more confusing then helpful in our case

    re_replace_in_file('\[external\]', '', './CGAL/html/annotated.html')

    # fix class/concept mismatch in generated pages
    relationship_pages=package_glob('./*/html/hasModels.html')
    relationship_pages.extend(package_glob('./*/html/generalizes.html'))
    relationship_pages.extend(package_glob('./*/html/refines.html'))
    for fn in relationship_pages:
        d = pq(filename=fn, parser='html')
        dts=d(".textblock .reflist dt")
        # no contents() on pyquery, do it the hard way
        # Note that in the following regular expression, the Struct did not appear in doxygen version 1.8.3
        # in hasModels.html, generalizes.html and refines.html, it is always Class. If this changes in
        # future versions of doxygen, the regular expression will be ready
        dts.each(lambda i: pq(this).html(re.sub("((Class )|(Struct ))", "Concept ", pq(this).html())))
        write_out_html(d, fn)

    # throw out nav-sync and the detailed description title
    all_pages=glob.glob('./*/html/*.html')
    for fn in all_pages:
        d = pq(filename=fn, parser='html')
        d('#nav-sync').hide()
        d('h2.groupheader').filter(lambda i: pq(this).text() == 'Detailed Description').remove()
        # TODO count figures
        write_out_html(d, fn)

    # remove %CGAL in navtree: this should be a fix in doxygen but for now it does not worth it
    re_replace_in_file('%CGAL','CGAL',glob.glob('./CGAL/html/navtree.js')[0])
    clean_doc()
    
    
if __name__ == "__main__":
    main()
