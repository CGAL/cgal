#!/usr/bin/env python2

import codecs
import re
import os
import sys

def main(argv):
    sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
    pattern = re.compile(r"\\package_listing{([^}]*)}")
    f = codecs.open(argv[1], 'r', encoding='utf-8')
    for line in f.readlines():
        match = pattern.match(line)
        if(match):
            pkg = match.group(1)
            index = pkg.find("/")
            if(index > 0):
                top_level = pkg[:index]
                lower_level = pkg[index+1:]
                filename="../" + top_level + "/doc/" + lower_level + "/PackageDescription.txt"
            else:
                filename="../" + pkg + "/doc/" + pkg + "/PackageDescription.txt"
            pkgdesc = codecs.open(filename, 'r', encoding='utf-8')
            do_print=False
            for l in pkgdesc.readlines():
                do_print = do_print or re.match(".*PkgDescriptionBegin.*", l)
                if(do_print):
                    sys.stdout.write(l)
                do_print = do_print and (not re.match(".*PkgDescriptionEnd.*", l))
        else:
            sys.stdout.write(line)

if __name__ == "__main__":
    main(sys.argv)
