#!/usr/bin/env python2

import codecs
import re
import os
import sys

def main(argv):
#    sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
    pattern = re.compile(r"\\package_listing{([^}]*)}")
    f = codecs.open(argv[1], 'r', encoding='utf-8')
    for line in f:
        match = pattern.match(line)
        if(match):
            pkg = match.group(1)
            index = pkg.find("/")
            if(index > 0):
                top_level = pkg[:index]
                lower_level = pkg[index+1:]
                filename="${CMAKE_SOURCE_DIR}/" + top_level + "/doc/" + lower_level + "/PackageDescription.txt"
            else:
                filename="${CMAKE_SOURCE_DIR}/" + pkg + "/doc/" + pkg + "/PackageDescription.txt"
            pkgdesc = codecs.open(filename, 'r', encoding='utf-8')
            do_print=False
            for l in pkgdesc:
                do_print = do_print or re.match(".*cgalPkgDescriptionBegin.*", l)
                if(do_print):
                    sys.stdout.write(l.encode('utf-8'))
                do_print = do_print and (not re.match(".*cgalPkgDescriptionEnd.*", l))
        else:
            sys.stdout.write(line.encode('utf-8'))

if __name__ == "__main__":
    main(sys.argv)
