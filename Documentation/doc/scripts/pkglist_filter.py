#!/usr/bin/env python

import codecs
import re
import os
import sys

SOURCE_DIR   = "${CGAL_ROOT}/"
BRANCH_BUILD = "${CGAL_BRANCH_BUILD}"

def make_doc_path(pkg, arg):
  if BRANCH_BUILD:
    return os.path.join(SOURCE_DIR, pkg, "doc", pkg, arg)
  else:
    return os.path.join(SOURCE_DIR, "doc", pkg, arg)

def main(argv):
#    sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
    pattern = re.compile(r"\\package_listing{([^}]*)}")
    f = codecs.open(argv[1], 'r', encoding='utf-8')
    for line in f:
        match = pattern.match(line)
        if(match):
            pkg = match.group(1)
            filename = make_doc_path(pkg, "PackageDescription.txt")
            pkgdesc = codecs.open(filename, 'r', encoding='utf-8')
            do_print=False
            for l in pkgdesc:
                do_print = do_print or re.match(".*cgalPkgDescriptionBegin.*", l)
                if(do_print):
                  if hasattr(sys.stdout, 'buffer'):
                    sys.stdout.buffer.write(l.encode('utf-8')) #python3
                  else:
                    sys.stdout.write(l.encode('utf-8')) #python2
                do_print = do_print and (not re.match(".*cgalPkgDescriptionEnd.*", l))
        else:
          if hasattr(sys.stdout, 'buffer'):
            sys.stdout.buffer.write(line.encode('utf-8')) #python3
          else:
            sys.stdout.write(line.encode('utf-8')) #python2

if __name__ == "__main__":
    main(sys.argv)
