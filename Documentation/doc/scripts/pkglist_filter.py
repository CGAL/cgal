#!/usr/bin/env python3
"""Filter package listings to emit package descriptions."""

import re
import os
import sys

SOURCE_DIR   = "${CGAL_ROOT}/"
BRANCH_BUILD = "${CGAL_BRANCH_BUILD}"

def make_doc_path(pkg, arg):
    """Return the path to a package documentation file."""
    if BRANCH_BUILD:
        return os.path.join(SOURCE_DIR, pkg, "doc", pkg, arg)
    return os.path.join(SOURCE_DIR, "doc", pkg, arg)

def main(argv):
    """Write filtered package descriptions to stdout in UTF-8."""
    pattern = re.compile(r"\\package_listing{([^}]*)}")
    with open(argv[1], 'r', encoding='utf-8') as infile:
        for line in infile:
            match = pattern.match(line)
            if match:
                pkg = match.group(1)
                filename = make_doc_path(pkg, "PackageDescription.txt")
                with open(filename, 'r', encoding='utf-8') as pkgdesc:
                    do_print = False
                    for pkg_line in pkgdesc:
                        do_print = do_print or re.match(
                            ".*cgalPkgDescriptionBegin.*", pkg_line
                        )
                        if do_print:
                            sys.stdout.buffer.write(pkg_line.encode('utf-8'))
                        do_print = do_print and (
                            not re.match(".*cgalPkgDescriptionEnd.*", pkg_line)
                        )
            else:
                sys.stdout.buffer.write(line.encode('utf-8'))

if __name__ == "__main__":
    main(sys.argv)
