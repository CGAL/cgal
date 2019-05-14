#!/usr/bin/env python

# This script creates the directory structure for a new package.
# Usage:
# cgal_create_package_dir.py Package_name [optional creation directory]

import sys
import os
import re
import argparse
import shutil

parser = argparse.ArgumentParser(
    description='Create directory structure for a new CGAL package.',
    epilog='A single directory named after the package, which contains ' +
           'the whole directory structure, is created at the creation path.')
parser.add_argument('packagename',
                    help='name of new CGAL package')
parser.add_argument('creationpath', nargs='?',
                    help='directory where package is created; ' +
                         'if omitted, the package directory is created ' +
                         'in the current directory')
args = parser.parse_args()
packagename = args.packagename
creationpath = args.creationpath

doxystring = \
r"""@INCLUDE = ${CGAL_DOC_PACKAGE_DEFAULTS}
PROJECT_NAME = "CGAL ${CGAL_DOC_VERSION} - Put title of project here"
INPUT        = ${CMAKE_SOURCE_DIR}/PACKAGENAME/doc/PACKAGENAME/ \
               ${CMAKE_SOURCE_DIR}/PACKAGENAME/include
"""

descrstring = \
r"""// PRETTY PACKAGE NAME should equal the project title in Doxyfile.in

/// \defgroup PkgPACKAGE PRETTY PACKAGE NAME Reference
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
\cgalPkgPicture{pkg-small.png}

\cgalPkgSummaryBegin
\cgalPkgAuthors{PACKAGE AUTHOR}
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

usermanstring = \
r"""namespace CGAL {
/*!

\mainpage User Manual
\anchor Chapter_PACKAGE_NAME
\anchor chaptermine
\cgalAutoToc
\author AUTHOR1, AUTHOR2

This chapter describes the ...

\section secmydefinitions Definitions

Section on definitions here ...

\section secmyexamples Examples

\subsection myFirstExample First Example

The following example shows ...

*/
} /* namespace CGAL */
"""

depsstring = \
r"""Manual
Kernel_23
STL_Extension
Algebraic_foundations
Circulator
Stream_support
"""

if re.match("^[A-Za-z_][A-Za-z0-9_]*$", packagename):

    if creationpath and (not creationpath == '.'):
        packagepath = os.path.join(creationpath, packagename)
    else:
        packagepath = packagename

    os.mkdir(packagepath)
    open(os.path.join(packagepath, 'dont_submit'), 'w').close()

    inclpath = os.path.join(packagepath, 'include', 'CGAL', packagename)
    os.makedirs(inclpath)

    os.mkdir(os.path.join(packagepath, 'src'))

    testpath = os.path.join(packagepath, 'test', packagename)
    os.makedirs(testpath)
    os.mkdir(os.path.join(testpath, 'data'))
    os.mkdir(os.path.join(testpath, 'include'))

    expath = os.path.join(packagepath, 'examples', packagename)
    os.makedirs(expath)
    os.mkdir(os.path.join(expath, 'data'))
    os.mkdir(os.path.join(expath, 'include'))
    open(os.path.join(expath, 'README'), 'w').close()

    demopath = os.path.join(packagepath, 'demo', packagename)
    os.makedirs(demopath)
    os.mkdir(os.path.join(demopath, 'data'))
    os.mkdir(os.path.join(demopath, 'include'))
    open(os.path.join(demopath, 'README'), 'w').close()

    benpath = os.path.join(packagepath, 'benchmark', packagename)
    os.makedirs(benpath)

    os.mkdir(os.path.join(packagepath, 'scripts'))
    os.mkdir(os.path.join(packagepath, 'developer_scripts'))

    infopath = os.path.join(packagepath, 'package_info', packagename)
    os.makedirs(infopath)
    open(os.path.join(infopath, 'copyright.txt'), 'w').close()
    open(os.path.join(infopath, 'description.txt'), 'w').close()
    open(os.path.join(infopath, 'license.txt'), 'w').close()
    open(os.path.join(infopath, 'long_description.txt'), 'w').close()
    open(os.path.join(infopath, 'maintainer'), 'w').close()

    docpath = os.path.join(packagepath, 'doc', packagename)
    os.makedirs(docpath)
    os.mkdir(os.path.join(docpath, 'CGAL'))
    os.mkdir(os.path.join(docpath, 'Concepts'))
    figpath = os.path.join(docpath, 'fig')
    os.mkdir(figpath)
    open(os.path.join(docpath, 'examples.txt'), 'w').close()

    usermanpath = os.path.join(docpath, (packagename + '.txt'))
    usermanfile = open(usermanpath, 'w')
    usermanfile.write(usermanstring)
    usermanfile.close()

    doxypath = os.path.join(docpath, ('Doxyfile.in'))
    doxyfile = open(doxypath, 'w')
    doxyfile.write(re.sub('PACKAGENAME', packagename, doxystring))
    doxyfile.close()

    descrpath = os.path.join(docpath, 'PackageDescription.txt')
    descrfile = open(descrpath, 'w')
    descrfile.write(descrstring)
    descrfile.close()

    depspath = os.path.join(docpath, ('dependencies'))
    depsfile = open(depspath, 'w')
    depsfile.write(depsstring)
    depsfile.close()

    # try to find figure pkg-small.png and copy it to figure path
    scriptdir = os.path.dirname(sys.argv[0])
    cgaldir = os.path.dirname(os.path.dirname(scriptdir))
    figfile = os.path.join(cgaldir, 'Documentation', 'doc',
        'Documentation', 'fig', 'pkg-small.png')

    if os.path.exists(figfile):
        shutil.copy(figfile, figpath)
else:
    sys.stderr.write("Error: Bad package name: " + packagename + '\n')
    sys.stderr.write("The package name should consist of:" + '\n')
    sys.stderr.write \
      ("letters, digits and underscores and not start with a digit." + '\n')
