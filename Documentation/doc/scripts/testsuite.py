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
import shutil
import sys
import subprocess
import os
import glob
import re
import operator
import datetime
from xml.dom.minidom import parseString
from pyquery import PyQuery as pq

def write_out_html(d, fn):
    with open(fn, 'w') as f:
        # this is the normal doxygen doctype, which is thrown away by pyquery
        f.write('<!DOCTYPE html>\n')
        f.write(d.html())

def count_errors_and_warnings(fn):
    with open(fn) as f:
        warn_count=0
        error_count=0
        for line in f:
            if re.match('^citelist.*[wW]arning.*$', line):
                continue
            if re.match('^.*[wW]arning.*$', line):
                warn_count=warn_count + 1
            if re.match('^.*[eE]rror.*$', line):
                error_count=error_count + 1
    return (warn_count, error_count)

def get_version():
    proc=subprocess.Popen(['git', 'rev-parse', 'HEAD'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    rev=proc.communicate()[0].strip()
    proc=subprocess.Popen(['git', 'log', '-n', '1', '--format=\"%ai\"', '--date=short'], stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True)
    date=proc.communicate()[0]
    #rev=subprocess.check_output(['git', 'rev-parse', 'HEAD'], universal_newlines=True)
    #date=subprocess.check_output(['git', 'log', '-n', '1', '--format=\"%ai\"', '--date=short'], universal_newlines=True)
    date=date[1:11]
    return (rev, date)

def write_report(args):
    page_header='''<html><head>
<title>CGAL Doxygen Manual Results</title>
<style type="text/css">
.test-results th {text-align: center;}
.test-results .warn-count {text-align: center;}
.test-results .error-count {text-align: center;}
.test-results .package-error {background-color: #FF8080;}
.test-results .package-good {background-color: #80FF80;}
.test-results .package-warnings {background-color: #FFFF80;}
p {margin-left:20px;}
body  {color: black; background-color: #C0C0D0; font-family: sans-serif;}
</style>
</head><body>
<h1 id="maintitle">Doxygen Manual Results</h1>'''
    page_footer='''<table class="test-results">
<tr>
<th>Package Name</th>
<th>Warnings</th>
<th>Errors</th>
</tr>
</table></body></html>'''
    
    if args.publish and args.do_copy_results:
      suffix=''
      if args.doxygen_version:
        suffix = ""+args.doxygen_version
      link="<a href=\"output/Manual/index.html\">Documentation built</a> with <a href=\"https://github.com/CGAL/doxygen\">our fork of Doxygen {_suffix}</a>\n".format(_suffix=suffix)
      suffix = ''
      if args.master_describe:
        suffix=args.master_describe
      link_master="\n<br><a href=\"master/Manual/index.html\">Documentation built</a> with <a href=\"https://github.com/doxygen/doxygen\">the master version of Doxygen {_suffix}</a> (buggy), so that we see progress/regression of doxygen development as far as CGAL is concerned.\n".format(_suffix=suffix)
      d = pq(page_header+link+"   "+link_master+page_footer)
    else:
      d = pq(page_header+page_footer)
    logs=sorted(glob.glob('./*.log'))
    err_war_sum=(0,0)
    for log in logs:
        res=count_errors_and_warnings(log)
        err_war_sum=tuple(map(operator.add, err_war_sum, res))
        status='class="package-error"'
        if res[0] == 0 and res[1] == 0:
            status='class="package-good"'
        elif res[0] != 0 and res[1] == 0:
            status='class="package-warnings"'
    
        basename=os.path.basename(log)
        pretty_name=basename[0:-4]
        new_row='''<tr {status}>
<td><a class="name" href="{basename}">{pretty_name}</a></td>
<td class="warn-count">{warn_count}
</td><td class="error-count">{err_count}</td></tr>'''.format(status=status, basename=basename, pretty_name=pretty_name, warn_count=str(res[0]), err_count=str(res[1]))

        d('.test-results').append(new_row)
    return (d, err_war_sum)

def main():
    parser = argparse.ArgumentParser(
    description='This script updates a checkout of cgal, purges the documentation, rebuilds it, creates an HTML summary of the resulting log files, and publishes the created files and logs.')

    parser.add_argument('--publish', metavar='/path/to/publish', help='Specify this argument if the results should be published.')
    parser.add_argument('--doc-log-dir', default='.', metavar='/path/to/cgal/build/dir/doc_log', help='The path of the documentation logs.')
    parser.add_argument('--master-dir', default='.', metavar='/path/to/cgal/build/master_dir/doc_output', help='The path to the master build documentation.')
    parser.add_argument('--output-dir', default='.', metavar='/path/to/cgal/build/dir/doc_output', help='The path to the build documentation')
    parser.add_argument('--diff', metavar='/path/to/diff', help='The path to the diff file.')
    parser.add_argument('--cgal-version', help='Path to a version.h file from the current release. If not specified use git hash instead.')
    parser.add_argument('--version-to-keep', help='indicates the number of release testsuites that should be kept at the publishing location.')
    parser.add_argument('--do-copy-results', action="store_true", help='Specify this argument if you want to copy the generated documentation into the publishing location.')
    parser.add_argument('--doxygen-version', default ='', help='Specify this argument if you want to add a version number to the name of the link to the documentation.')
    parser.add_argument('--master-describe', default ='', help='Specify this argument if you want to add a suffix to the name of the link to the doxygen master documentation.')
    
    args = parser.parse_args()
    
    os.chdir(args.doc_log_dir)

    d, sum=write_report(args)
    if args.cgal_version:
      version_string="CGAL-"+args.cgal_version
      version_date=datetime.datetime.now().strftime("%Y-%m-%d")
    else:
      version_string,version_date=get_version()
      version_string="Revision-"+version_string

    title=d('#maintitle')
    title.text(title.text() + ' for ' + version_string)
    write_out_html(d, './index.html')

    # does the diff exist ?
    diff='n/a'
    if args.diff:
        diff_file=args.diff
        if not os.path.isfile(diff_file):
            sys.stderr.write('Diff file ' + diff_file + ' is not a file. Cannot diff.\n')
            sys.exit(1)    
    if args.publish:
        if args.publish.endswith('/'):
            publish_dir=args.publish
        else:
            publish_dir=args.publish + '/'

        if not os.path.isdir(publish_dir):
            sys.stderr.write('Publish dir ' + publish_dir + ' is not a directory. Cannot publish.\n')
            sys.exit(1)

        if os.path.isdir(publish_dir + 'log' + version_string):
            sys.stderr.write('Logs for this revision have already been publish under: ' 
                             + publish_dir + 'log' + version_string + ' Cannot publish.\n')
            sys.exit(1)
        log_target=publish_dir + version_string
        # does the index file exist? if not write out a skeleton
        try:
            with open(publish_dir + 'index.html') as f: pass
        except IOError as e:
            print('No index.html in the publish directory found. Writing a skeleton.')               
            with open(publish_dir + 'index.html', 'w') as f:
                f.write('''<!DOCTYPE html>
<style type="text/css">
.rev-table th {text-align: center;}
.rev-table td {text-align: center;}
table {margin-left:40px;}
body  {color: black; background-color: #C0C0D0; font-family: sans-serif;}
</style>
<html><head><title>Manual Testsuite Overview</title></head>
<body><h1>Overviewpage of the Doxygen Manual Testsuite</h1>
<table  border="1" cellspacing="2" cellpadding="5" id="revisions" class="rev-table">
<tr><th>Revision</th><th>Date</th><th>Warnings</th><th>Errors</th><th>Diff with doxygen master</th></tr></table></body></html>''')

        with open(diff_file, 'r') as myfile:
          diff=myfile.read()
        if not diff:
          diff='none'
        else:
          diff='<a href="{log_path}/diff.txt">Diff between {test_version} and {master_version}.</a>'.format(
          log_path=version_string, test_version=args.doxygen_version, master_version=args.master_describe)
        d=pq(filename=publish_dir + 'index.html',parser="html")
        revs=d('#revisions tr')
        new_row='<tr><td><a href="{revision}/index.html">{revision}</a></td><td>{date}</td><td>{warnings}</td><td>{errors}</td><td>{diffs}</td></tr>'.format(
            revision=version_string, date=version_date, warnings=sum[0], errors=sum[1], diffs=diff)
        revs.eq(0).after(new_row)
        if args.version_to_keep:
          nb_items=len(revs)
          for k in range(int(args.version_to_keep),nb_items):
            dir_to_remove=revs.eq(k).text().split()[0]
            if os.access(publish_dir + dir_to_remove, os.W_OK):
                shutil.rmtree(publish_dir + dir_to_remove)
            else:
                sys.stderr.write("Warning: the directory " + publish_dir + dir_to_remove + " does not exist or is not writable!\n")
            revs.eq(k).remove()
        write_out_html(d, publish_dir + 'index.html')
        try:
          #copy log files
          shutil.copytree('.', log_target)
          #copy diff
          shutil.copyfile(diff_file, log_target+'/diff.txt')
          try:
            #copy documentation
            if args.do_copy_results:
              tgt=os.path.join(log_target, 'output')
              shutil.copytree(args.output_dir, tgt, symlinks=True)
              os.symlink("../MathJax", os.path.join(log_target, 'MathJax'))
          except:
            sys.stderr.write("Error while copying documentation\n")
            raise
          if args.master_dir:
            try:
              #copy documentation from master
              if args.do_copy_results:
                tgt=os.path.join(log_target, 'master')
                shutil.copytree(args.master_dir, tgt, symlinks=True)
            except:
              sys.stderr.write("Error while copying master documentation\n")
              raise

        except:
          sys.stderr.write("Error while writing to "+log_target+". Does it already exists?\n")
          raise
        
if __name__ == "__main__":
    main()
