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

xmlstring = \
"""<project>
 <name>PROJECT NAME</name>
 <input>../PACKAGENAME/doc/PACKAGENAME</input>
 <doxygen>
   <string name="STRIP_FROM_PATH">../PACKAGENAME/doc/PACKAGENAME</string>
   <string name="STRIP_FROM_INC_PATH">../PACKAGENAME/doc/PACKAGENAME/</string>
   <string name="GENERATE_TAGFILE">./tags/PACKAGENAME.tag</string>
   <string name="IMAGE_PATH">../PACKAGENAME/doc/PACKAGENAME/fig</string>
   <string name="EXAMPLE_PATH">../PACKAGENAME/examples</string>
   <list name="PACKAGE_REFERENCES">
       <item>PACKAGENAME</item>
   </list>
 </doxygen>
</project>
"""

descrstring = \
r"""/// \defgroup PkgPACKAGE PACKAGE NAME Reference
/// \defgroup PkgPACKAGEConcepts Concepts
/// \ingroup PkgPACKAGE

/// \defgroup PkgPACKAGEAlgorithmClasses Algorithm Classes
/// \ingroup PkgPACKAGE

/// \defgroup PkgPACKAGETraitsClasses Traits Classes
/// \ingroup PkgPACKAGE

/// \defgroup PkgPACKAGEMiscellaneous Miscellaneous
/// \ingroup PkgPACKAGE

/*!
\addtogroup PkgPACKAGE
\todo check generated documentation

\cgalPkgDescriptionBegin{PACKAGE NAME,PkgPACKAGESummary}
\cgalPkgPicture{cdt2d-small.png}

\cgalPkgSummaryBegin
\cgalPkgAuthor{PACKAGE AUTHOR}
\cgalPkgDesc{PACKAGE DESCRIPTION.
The package provides ... }
\cgalPkgManuals{Chapter_PACKAGE_NAME,PkgPACKAGE}
\cgalPkgSummaryEnd

\cgalPkgShortInfoBegin
\cgalPkgSince{X.X}
\cgalPkgDependsOn{\ref PkgDEPENDENCY}
\cgalPkgBib{cgal:x-x}
\cgalPkgLicense{\ref licensesGPL "GPL"}
\cgalPkgDemo{DEMO 1,demo1.zip}
\cgalPkgDemo{DEMO 2,demo2.zip}
\cgalPkgShortInfoEnd

\cgalPkgDescriptionEnd

*/
"""

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

    infopath = os.path.join(packagename, 'package_info', packagename)
    os.makedirs(infopath)
    open(os.path.join(infopath, 'copyright.txt'), 'w').close()
    open(os.path.join(infopath, 'description.txt'), 'w').close()
    open(os.path.join(infopath, 'license.txt'), 'w').close()
    open(os.path.join(infopath, 'long_description.txt'), 'w').close()
    open(os.path.join(infopath, 'maintainer'), 'w').close()

    docpath = os.path.join(packagename, 'doc', packagename)
    os.makedirs(docpath)
    os.mkdir(os.path.join(docpath, 'CGAL'))
    os.mkdir(os.path.join(docpath, 'Concepts'))
    os.mkdir(os.path.join(docpath, 'fig'))
    open(os.path.join(docpath, 'examples.txt'), 'w').close()
    open(os.path.join(docpath, (packagename + '.txt')), 'w').close()

    xmlpath = os.path.join(docpath, (packagename + '.xml'))
    xmlfile = open(xmlpath, 'w')
    xmlfile.write(re.sub('PACKAGENAME', packagename, xmlstring))
    xmlfile.close()

    descrpath = os.path.join(docpath, 'PackageDescription.txt')
    descrfile = open(descrpath, 'w')
    descrfile.write(descrstring)
    descrfile.close()

else:
    print "Error: Bad package name:", packagename
    print "The package name should consist of:"
    print "letters, digits and underscores and not start with a digit."
