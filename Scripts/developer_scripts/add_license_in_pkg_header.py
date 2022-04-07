#!/usr/bin/env python

# Insert in all the header file of a package the include directive
# "#include <CGAL/license/${package_name}.h" if it is not already there
# There are two arguments to the script as some packages are split into
# subdirectories (Algebraic_kernel_for_circles, ...)

from sys import argv
import os
from os.path import join
import codecs
import re

def add_license_include_in_file(package_name, fname):
  # first, see if the include directive is already there
  with codecs.open(fname, encoding='utf-8') as f:
    if any(re.search("#include\s*<CGAL/license/", line) for line in f):
      return # include directive already there

  #match only file under a GPL license
  with codecs.open(fname, encoding='utf-8') as f:
    if not any(re.search("SPDX-License-Identifier:.*[ (]GPL", line) for line in f):
      return # include directive already there


  # include directive not already there
  inserted = False
  with codecs.open(fname, encoding='utf-8') as f:
    out_fname = fname + ".tmp"
    out = codecs.open(out_fname, "w", encoding='utf-8')
    for line in f:
      if not inserted and re.search("#\s*define\s+CGAL_.*_H\s*$", line):
        out.write(line+"\n")
        out.write("#include <CGAL/license/"+package_name+".h>\n\n")
        inserted=True
      else:
        out.write(line)
    out.close()
    f.close()
    os.remove(fname)
    os.rename(out_fname, fname)
    if not inserted:
      print("Warning: file "+fname+" was not modified (no CGAL_*_H defined)")

if len(argv)==1:
  print("Usage: "+argv[0]+" Package_directory [Package_name=Package_directory]\n")
else:
  package_dir=argv[1]
  if len(argv)==3:
    package_name=argv[2]
  else:
    package_name=package_dir

  for root, dirs, files in os.walk(argv[1]+'/include/CGAL'):
    for f in files:
      if f.endswith('.h'):
        add_license_include_in_file(package_name,os.path.join(root,f))

