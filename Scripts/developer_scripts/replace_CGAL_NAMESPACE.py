#!/usr/bin/python

import os
import fileinput
from sys import argv
from sys import exit
from sys import stderr

if len(argv)!=2:
  stderr.write("Usage: python "+argv[0]+" directory_name.\n")
  stderr.write("Action:\nReplace CGAL_BEGIN_NAMESPACE and CGAL_END_NAMESPACE by \'namespace CGAL {\'\n")
  stderr.write("and \'} //namespace CGAL\' respectively in all .h and .cpp files in the directory \n")
  stderr.write("given as parameter.\n")

for dirname, dirnames, filenames in os.walk(argv[1]):
    for filename in filenames:
        l=len(filename)
        if ( l>2 and filename[(l-2):l]==".h") or ( l> 4 and filename[(l-4):l]==".cpp"):
          fname=os.path.join(dirname, filename)
          for lines in fileinput.FileInput(fname, inplace=1):
            lines = lines.replace("CGAL_BEGIN_NAMESPACE","namespace CGAL {")
            lines = lines.replace("CGAL_END_NAMESPACE","} // end namespace CGAL")
            print lines,
