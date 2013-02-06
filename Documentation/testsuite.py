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
            if re.match('^citelist.*warning.*$', line):
              continue
            if re.match('^.*: warning: .*$', line):
                warn_count=warn_count + 1
            if re.match('^.*: error: .*$', line):
                error_count=error_count + 1
    return (warn_count, error_count)

def update():
    subprocess.call(['git', 'pull'])

def integration_test():
  subprocess.call(['git', 'fetch'])
  subprocess.call(['git', 'checkout','origin/integration'])
  subprocess.call(['git', 'reset','--hard','origin/integration'])

def purge_doc():
    for log in glob.glob('./log/*.*'):
        os.remove(log)

    for tag in glob.glob('./tags/*.*'):
        os.remove(tag)

    for output in glob.glob('./output/CGAL.CGAL*'):
        shutil.rmtree(output)

def run_doxyassist(doxyassist, doxygen):
    subprocess.call([sys.executable,doxyassist, '--debug', '--doxygen', doxygen, 'doxyassist.xml'])
    
def get_version():
    rev=subprocess.check_output(['git', 'rev-parse', 'HEAD'], universal_newlines=True)
    date=subprocess.check_output(['git', 'log', '-n', '1', '--format=\"%ai\"', '--date=short'], universal_newlines=True)
    date=date[1:11]
    return (rev, date)

def get_cgal_version(path_to_version_h):
    f=file(path_to_version_h)
    for line in f.readlines():
        m = re.match('^#define CGAL_VERSION (.*)',line)
        if m:
          return "CGAL-"+m.group(1)
    return "Unknown Version"


def write_report():
    d = pq('''<html><head>
<title>CGAL Doxygen Manual Results</title>
<style type="text/css">
.test-results th {text-align: center;}
.test-results .warn-count {text-align: center;}
.test-results .error-count {text-align: center;}
.test-results .package-error {background-color: #FF8080;}
.test-results .package-good {background-color: #80FF80;}
.test-results .package-warnings {background-color: #FFFF80;}
p {margin-left:20px;}
body {background-image:url("images/back40.gif");}
</style>
</head><body>
<h1 id="maintitle">Doxygen Manual Results</h1><table class="test-results">
<tr>
<th>Package Name</th>
<th>Warnings</th>
<th>Errors</th>
</tr>
</table></body></html>''')
    logs=sorted(glob.glob('./log/CGAL.CGAL*-error.log'))
    err_war_sum=(0,0)
    for log in logs:
        res=count_errors_and_warnings(log)
        err_war_sum=tuple(map(operator.add, err_war_sum, res))
        status='class="package-errors"'
        if res[0] == 0 and res[1] == 0:
            status='class="package-good"'
        elif res[0] != 0 and res[1] == 0:
            status='class="package-warnings"'
    
        basename=os.path.basename(log)
        pretty_name=basename[5:-10]
        new_row='''<tr {status}>
<td><a class="name" href="{basename}">{pretty_name}</a></td>
<td class="warn-count">{warn_count}
</td><td class="error-count">{err_count}</td></tr>'''.format(status=status, basename=basename, pretty_name=pretty_name, warn_count=str(res[0]), err_count=str(res[1]))

        d('.test-results').append(new_row)
    return (d, err_war_sum)

def main():
    parser = argparse.ArgumentParser(
    description='This script updates a checkout of cgal, purges the documentation, rebuilds it, creates an HTML summary of the resulting log files, and publishes the created files and logs.')
    parser.add_argument('--doxyassist', default='/usr/bin/doxyassist.py', metavar='/path/to/doxyassist.py')
    parser.add_argument('--doxygen', default='/usr/bin/doxygen', metavar='/path/to/doxygenbinary', help='the doxygen binary', )
    parser.add_argument('--documentation', default='.', metavar='/path/to/cgal/Documentation', help='The path to the Documentation dir of the git checkout you would like to test.')
    parser.add_argument('--publish', metavar='/path/to/publish', help='Specify this argument if the results should be published.')
    parser.add_argument('--do-update', action="store_true", help='Specify this argument if you want to do a version control update.')
    parser.add_argument('--test-integration', action="store_true", help='Specify this argument if you want to switch to integration and use the latest version.')
    parser.add_argument('--do-purge-rebuild', action="store_true", help='Specify this argument if you want to actually rebuild the documentation. Just write the report if not specified.')
    parser.add_argument('--cgal-version', help='Path to a version.h file from the current release. If not specified use git hash instead.')
    
    args = parser.parse_args()
    
    os.chdir(args.documentation)
    if args.do_update:
        update()

    if args.test_integration:
      integration_test()

    if args.do_purge_rebuild:
        doxyassist="".join(args.doxyassist)
        doxygen="".join(args.doxygen)
        purge_doc()
        # two runs are required, one to build the tags, the next to actually use them
        run_doxyassist(doxyassist, doxygen)
        run_doxyassist(doxyassist, doxygen)
        subprocess.call([sys.executable,'./html_output_post_processing.py', '--output', './output'])

    d, sum=write_report()
    version_string,version_date=get_version()
    if args.cgal_version:
      version_string=get_cgal_version(args.cgal_version)
    else:
      version_string="Revision "+version_string

    title=d('#maintitle')
    title.text(title.text() + ' for ' + version_string)
    write_out_html(d, './log/index.html')
    
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
body {background-image:url("images/back40.gif");}
</style>
<html><head><title>Manual Testsuite Overview</title></head>
<body><h1>Overviewpage of the Doxygen Manual Testsuite</h1><table id="revisions" class="rev-table"><tr><th>Revision</th><th>Date</th><th>Warnings</th><th>Errors</th></tr></table></body></html>''')

        d=pq(filename=publish_dir + 'index.html',parser="html")
        revs=d('#revisions tr')
        new_row='<tr><td><a href="log{revision}/testsuite.html">{revision}</a></td><td>{date}</td><td>{warnings}</td><td>{errors}</td></tr>'.format(
            revision=version_string, date=version_date, warnings=sum[0], errors=sum[1])
        revs.eq(0).after(new_row)
        write_out_html(d, publish_dir + 'index.html')
        shutil.copytree('./log', publish_dir + 'log-' + version_string)
        
if __name__ == "__main__":
    main()
