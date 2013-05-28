#!/usr/bin/python

# This script creates the directory structure for a new package.
# Usage: cgal_create_package_dir.py Package_name

import sys
import os
import re
import argparse

parser = argparse.ArgumentParser(
    description='Create directory structure for a new CGAL package.')
parser.add_argument('packagename',
                    help='name of new CGAL package')
args = parser.parse_args()
packagename = args.packagename

if re.match("^[A-Za-z_][A-Za-z0-9_]*$", packagename):
    os.mkdir(packagename)
    open(os.path.join(packagename, 'dont_submit'), 'w').close()

    os.makedirs(os.path.join(packagename, 'include', 'CGAL'))

    os.mkdir(os.path.join(packagename, 'src'))

    testpath = os.path.join(packagename, 'test', packagename)
    os.makedirs(testpath)
    os.mkdir(os.path.join(testpath, 'data'))
    os.mkdir(os.path.join(testpath, 'include'))

    expath = os.path.join(packagename, 'examples', packagename)
    os.makedirs(expath)
    os.mkdir(os.path.join(expath, 'data'))
    os.mkdir(os.path.join(expath, 'include'))
    open(os.path.join(expath, 'README'), 'w').close()

    demopath = os.path.join(packagename, 'demo', packagename)
    os.makedirs(demopath)
    os.mkdir(os.path.join(demopath, 'data'))
    os.mkdir(os.path.join(demopath, 'include'))
    open(os.path.join(demopath, 'README'), 'w').close()

    benpath = os.path.join(packagename, 'benchmark', packagename)
    os.makedirs(benpath)

    os.mkdir(os.path.join(packagename, 'scripts'))
    os.mkdir(os.path.join(packagename, 'developer_scripts'))

    infopath = os.path.join(packagename,'package_info', packagename)
    os.makedirs(infopath)
    open(os.path.join(infopath, 'copyright.txt'), 'w').close()
    open(os.path.join(infopath, 'description.txt'), 'w').close()
    open(os.path.join(infopath, 'license.txt'), 'w').close()
    open(os.path.join(infopath, 'long_description.txt'), 'w').close()
    open(os.path.join(infopath, 'maintainer'), 'w').close()

else:
    print "Error: Bad package name:", packagename
    print "The package name should consist of:"
    print "letters, digits and underscores and not start with a digit."
