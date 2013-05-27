#!/usr/bin/python

# This script creates the directory structure for a new package.
# Usage: cgal_create_package_dir.py Package_name

import sys
import os
import re

if len(sys.argv) != 2:
    print 'Usage:', sys.argv[0], 'Package_name'
elif re.match("^[A-Za-z_][A-Za-z0-9_]*$", sys.argv[1]):
    packagename = sys.argv[1]
    #print 'Creating package directory named', packagename
    os.mkdir(packagename)
    open(packagename + '/dont_submit', 'w').close()

    os.makedirs(packagename + '/include/CGAL')

    os.mkdir(packagename + '/src')

    testdirname = packagename + '/test/' + packagename
    os.makedirs(testdirname)
    os.mkdir(testdirname + '/data')
    os.mkdir(testdirname + '/include')

    exdirname = packagename + '/examples/' + packagename
    os.makedirs(exdirname)
    os.mkdir(exdirname + '/data')
    os.mkdir(exdirname + '/include')
    open(exdirname + '/README', 'w').close()

    demodirname = packagename + '/demo/' + packagename
    os.makedirs(demodirname)
    os.mkdir(demodirname + '/data')
    os.mkdir(demodirname + '/include')
    open(demodirname + '/README', 'w').close()

    bendirname = packagename + '/benchmark/' + packagename
    os.makedirs(bendirname)

    os.mkdir(packagename + '/scripts')
    os.mkdir(packagename + '/developer_scripts')
  
    infodirname = packagename + '/package_info/' + packagename
    os.makedirs(infodirname)
    open(infodirname + '/copyright.txt', 'w').close()
    open(infodirname + '/description.txt', 'w').close()
    open(infodirname + '/license.txt', 'w').close()
    open(infodirname + '/long_description.txt', 'w').close()
    open(infodirname + '/maintainer', 'w').close()

else:
    print "Error: Bad package name"
    print "Error: It should consist of letters, digits", \
          "and underscores and not start with a digit"
