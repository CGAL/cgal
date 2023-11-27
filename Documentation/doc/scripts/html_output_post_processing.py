#!/usr/bin/env python
# Copyright (c) 2012 GeometryFactory (France). All rights reserved.
# All rights reserved.
#
# This file is part of CGAL (www.cgal.org).
#
# $URL$
# $Id$
# SPDX-License-Identifier: GPL-3.0-or-later
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
    title = d("title")
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
    f.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "https://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n')
    f.write('<html xmlns=\"http://www.w3.org/1999/xhtml\">')
    if d.html() is not None:
      f.write(d.html(method='html'))
    f.write('\n')
    f.write('</html>\n')
    f.close()

def package_glob(target):
    return [x for x in glob.glob(target) if not os.path.join(os.path.join('.','Manual'),'')  in x]

# remove duplicate files
def clean_doc():
    duplicate_files=list(package_glob('./*/jquery.js'))
    duplicate_files.extend(package_glob('./*/dynsections.js'))
    duplicate_files.extend(package_glob('./*/resize.js'))
    duplicate_files.extend(package_glob('./*/cgal_stylesheet.css'))
    # kill _all_, including the one in CGAL tabs.css files
    duplicate_files.extend(glob.glob('./*/tabs.css'))
    # left-over by doxygen?
    duplicate_files.extend(package_glob('./*/bib2xhtml.pl'))
    duplicate_files.extend(package_glob('./*/cgal_manual.bib'))
    duplicate_files.extend(package_glob('./*/citelist.doc'))
    duplicate_files.extend(package_glob('./*/doxygen.bst'))
    duplicate_files.extend(package_glob('./*/geom.bib'))

    for fn in duplicate_files:
        os.remove(fn)

# from https://stackoverflow.com/a/1597755/105672
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
        f.close()
        os.remove(fname)
        os.rename(out_fname, fname)

def re_replace_first_in_file(pat, s_after, fname):
    # first, see if the pattern is even in the file.
    with codecs.open(fname, encoding='utf-8') as f:
        if not any(re.search(pat, line) for line in f):
            return # pattern does not occur in file so we are done.

    found=False
    # pattern is in the file, so perform replace operation.
    with codecs.open(fname, encoding='utf-8') as f:
        out_fname = fname + ".tmp"
        out = codecs.open(out_fname, "w", encoding='utf-8')
        for line in f:
            if not found and re.search(pat, line):
              out.write(re.sub(pat, s_after, line))
              found=True
            else:
              out.write(line)
        out.close()
        f.close()
        os.remove(fname)
        os.rename(out_fname, fname)

def is_concept_file(filename):
  if not path.exists(filename):
    return False;
  file_content = codecs.open(filename, 'r', encoding='utf-8')
  d = pq(file_content.read(),parser="html")
  ident = d('#CGALConcept')
  return ident.size() == 1

def rearrange_icon(i, dir_name):
    icon = pq(this)
    if icon.attr("class") == "icon-class":
        parser=pq(this).parent().parent()
        for link_class in ['a.el', 'a.elRef']:
            links=parser(link_class)
            if links.size()>0 and is_concept_file(path.join(dir_name, pq(links[0]).attr("href"))):
                icon.attr("class","icon-concept")

###############################################################################
############################## Figure Numbering ###############################
###############################################################################

## A helper for collecting the figure anchors in a package
class figure_anchor_info:
  def __init__(self,pkg_id,global_anchor_map):
    self.next_index=1
    self.global_anchor_map=global_anchor_map
    self.pkg_id=str(pkg_id)

## Collects the figure anchors in a package
def collect_figure_anchors(i,infos):
  anchor_name=pq(this).attr('id')
  if re.match("fig__.+",anchor_name) != None:
    infos.global_anchor_map[anchor_name]=infos.pkg_id+"."+str(infos.next_index)
    infos.next_index+=1

## Changes an anchor name using the name in global_anchor_map
def update_figure_ref(i,global_anchor_map):
  link = pq(this)
  link_name=link.text()
  if re.match("fig__.+",link_name) != None:
    #replace the link only if it was collected
    if link_name in global_anchor_map:
      link.text( "Figure "+str(global_anchor_map[link_name]) )
    else:
      stderr.write("Error: Figure numbering; "+link_name+" has not been collected\n")

## number figures and reference to figures
def automagically_number_figures():
  #collect the list of packages in the package overview page,
  #respecting the order of that page
  all_packages=[]
  if not path.isfile("./Manual/packages.html"):
    stderr.write("Error: Figure numbering; ./Manual/packages.html does not exist\n")
    return


  file_content = codecs.open("./Manual/packages.html", 'r', encoding='utf-8')
  d = pq(file_content.read(),parser="html")
  for el in d('a.elRef'):
    text = pq(el).attr('href')
    if text.find("index.html")!=-1:
      re_pkg_index=re.compile("\.\./([A-Z_a-z0-9]+)/index\.html")
      res=re_pkg_index.match(text)
      if res:
        all_packages.append(res.group(1))
      else:
        stderr.write("Error: Figure numbering; skipping "+text+"\n")
  #collect all the figure anchors and associate them a unique id
  pkg_id=1 # the id of a package
  global_anchor_map={} #map a figure anchor to it unique name: pkg_id+.+fig_nb
  for pkg_name in all_packages:
    # collect first in the user manual
    userman=os.path.join(pkg_name,"index.html")
    all_pkg_files=glob.glob(os.path.join(pkg_name,"*.html"))
    all_pkg_files.remove(userman)
    for fname in [userman]+all_pkg_files:
      infos=figure_anchor_info(pkg_id, global_anchor_map)
      file_content = codecs.open(fname, 'r', encoding='utf-8')
      d = pq(file_content.read(), parser="html")
      d('a.anchor').each( lambda i: collect_figure_anchors(i,infos) )
    pkg_id+=1

  #Figure link dev Manual
  for fname in glob.glob("Manual/*.html"):
    infos=figure_anchor_info(0, global_anchor_map)
    file_content = codecs.open(fname, 'r', encoding='utf-8')
    d = pq(file_content.read(),parser="html")
    d('a.anchor').each( lambda i: collect_figure_anchors(i,infos) )

  #replace each link to a figure by its unique id
  all_files=glob.glob("*/*.html")
  for fname in all_files:
    #consider only a html files containing a figure ref.
    with codecs.open(fname, encoding='utf-8') as f:
        if not any(re.search("fig__", line) for line in f):
            continue # pattern does not occur in file so we are done.
    file_content = codecs.open(fname, 'r', encoding='utf-8')
    d = pq(file_content.read(), parser="html")
    d('a.el').each( lambda i: update_figure_ref(i,global_anchor_map) )
    d('a.elRef').each( lambda i: update_figure_ref(i,global_anchor_map) )
    file_content.close()
    write_out_html(d, fname)

###############################################################################
###############################################################################
###############################################################################

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

    #workaround CGAL link on the main page
    re_replace_in_file("<a class=\"elRef\" href=\"../.+/namespaceCGAL.+html\">CGAL</a>", "CGAL", "./Manual/index.html")

    #workaround issue with operator<< in pyquery
    all_pages=glob.glob('*/*.html')
    for f in all_pages:
      re_replace_in_file("operator<<\(\)", "operator&lt;&lt;()", f)

    # number figure
    automagically_number_figures()

    #replace icons with CGAL colored ones
    annotated_files=package_glob('./*/annotated.html')
    for fn in annotated_files:
      re_replace_in_file("<span class=\"icon\">N</span>", "<span class=\"icon-namespace\">N</span>", fn)
      re_replace_in_file("<span class=\"icon\">C</span>", "<span class=\"icon-class\">C</span>", fn)
      dir_name=path.dirname(fn)
      file_content = codecs.open(fn, 'r', encoding='utf-8')
      d = pq(file_content.read(), parser="html")
      span_tags = d('table.directory tr span')
      span_tags.each(lambda i: rearrange_icon(i, dir_name))
      file_content.close()
      write_out_html(d,fn)
    class_files=list(package_glob('./*/class*.html'))
    class_files.extend(package_glob('./*/struct*.html'))
    for fn in class_files:
        file_content = codecs.open(fn, 'r', encoding='utf-8')
        d = pq(file_content.read(), parser="html")
        ident = d('#CGALConcept')
        if ident.size() == 1:
            conceptify(d);
            # in case of a second pass don't process this again
            d.remove("#CGALConcept")
        # there is a doxygen bug that prevents the correct linkage of the CGAL breadcrumb
        ident = d('#nav-path .navelem').eq(0).children().eq(0)
        if ident and ident.attr('href') == 'namespaceCGAL.html':
            ident.attr('href', '../Manual/namespaceCGAL.html')
        file_content.close()
        write_out_html(d, fn)

    namespace_files=package_glob('./*/namespace*.html')
    for fn in namespace_files:
        file_content = codecs.open(fn, 'r', encoding='utf-8')
        d = pq(file_content.read(), parser="html")
        ident = d('#CGALConceptNS')
        if ident.size() == 1:
            conceptify_ns(d);
            d.remove("#CGALConceptNS")
        file_content.close()
        write_out_html(d, fn)

    # in a group we only need to change the nested-classes
    group_files=package_glob('./*/group*Concepts*.html')
    for fn in group_files:
        file_content = codecs.open(fn, 'r', encoding='utf-8')
        d = pq(file_content.read(), parser="html")
        conceptify_nested_classes(d)
        file_content.close()
        write_out_html(d, fn)

    # fix up Files
    files_files=package_glob('./*/files.html')
    for fn in files_files:
        file_content = codecs.open(fn, 'r', encoding='utf-8')
        d = pq(file_content.read(), parser="html")
        table = d("table.directory")
        row_id=table("td.entry").filter(lambda i: pq(this).text() == 'Concepts').parent().attr('id')
        if row_id != None:
            # figure out the rowid and then drop everything from the table that matches
            table("tr").filter(lambda i: re.match(row_id + '*', pq(this).attr('id'))).remove()
            file_content.close()
            write_out_html(d, fn)

    #Rewrite the code for index trees images

    filesjs_files=package_glob('./*/files.js')
    for fn in filesjs_files:
      re_replace_in_file('^.*\[ "Concepts",.*$', '', fn)

    #Rewrite the path of some images
    re_replace_in_file("'src','ftv2",
                       "'src','../Manual/ftv2",
                       os.path.join('Manual','dynsections.js') )

    # external is placed by doxygen to mark a class from a tagfile, this
    # is more confusing then helpful in our case
    if path.isfile(os.path.join('Manual','annotated.html')):
      re_replace_in_file('\[external\]', '', os.path.join('Manual','annotated.html'))
    else:
      stderr.write("Error: ./Manual/annotated.html does not exists\n")
    # fix class/concept mismatch in generated pages
    relationship_pages=list(package_glob('./*/hasModels.html'))
    relationship_pages.extend(package_glob('./*/generalizes.html'))
    relationship_pages.extend(package_glob('./*/refines.html'))
    for fn in relationship_pages:
        file_content = codecs.open(fn, 'r', encoding='utf-8')
        d = pq(file_content.read(), parser="html")
        dts=d(".textblock .reflist dt")
        # no contents() on pyquery, do it the hard way
        # Note that in the following regular expression, the Struct did not appear in doxygen version 1.8.3
        # in hasModels.html, generalizes.html and refines.html, it is always Class. If this changes in
        # future versions of doxygen, the regular expression will be ready
        dts.each(lambda i: pq(this).html(re.sub("((Class )|(Struct ))", "Concept ", pq(this).html())))
        file_content.close()
        write_out_html(d, fn)

    # throw out nav-sync
    all_pages=glob.glob('./*/*.html')
    for fn in all_pages:
        file_content = codecs.open(fn, 'r', encoding='utf-8')
        d = pq(file_content.read(), parser="html")
        d('#nav-sync').hide()
        # TODO count figures
        file_content.close()
        write_out_html(d, fn)

    # remove %CGAL in navtree: this should be a fix in doxygen but for now it does not worth it
    re_replace_first_in_file('%CGAL','CGAL',glob.glob('./Manual/navtree.js')[0])
    clean_doc()

    #remove links to CGAL in the bibliography
    citelist_files=package_glob('./*/citelist.html')
    for fn in citelist_files:
        re_replace_in_file('<a class=\"el\" href=\"namespaceCGAL.html\">CGAL</a>', 'CGAL', fn)

    #add a section for Inherits from
    class_and_struct_files=list(package_glob('./*/class*.html'))
    class_and_struct_files.extend(package_glob('./*/struct*.html'))
    for fn in class_and_struct_files:
        re_replace_first_in_file(r'<p>Inherits\s*(.*)</p>', r'<h2 class="groupheader">Inherits from</h2><p>\1</p>', fn)

    # remove class name in Definition section if there is no default template
    # parameter documented
    for fn in class_and_struct_files:
        file_content = codecs.open(fn, 'r', encoding='utf-8')
        d = pq(file_content.read(), parser="html")
        for el in d('h3'):
          text = pq(el).text()
          if text[0:9]=="template<" and text.find('=')==-1:
            pq(el).remove()
        file_content.close()
        write_out_html(d, fn)

    #add a canonical link to all pages
    all_pages=glob.glob('*/*.html')
    for f in all_pages:
      url_f=os.path.split(f)
      url_f=url_f[0]+"/"+url_f[1]
      canonical_link="<link rel=\"canonical\" href=\"https://doc.cgal.org/latest/"+url_f+"\"/>\n"
      re_replace_first_in_file(r'<head>', r'<head>\n'+canonical_link, f)
    ## special case for how_to_cite.html
    canonical_link="<link rel=\"canonical\" href=\"https://doc.cgal.org/latest/Manual/how_to_cite.html\"/>\n"
    re_replace_first_in_file(r'<body>', r'<head>\n'+canonical_link+"</head>\n<body>", os.path.join("Manual","how_to_cite.html"))

if __name__ == "__main__":
    main()
