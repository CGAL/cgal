// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-3 $
// release_date  : $CGAL_Date: 1999/03/08 $
//
// file          : config/testfiles/CGAL_CFG_NO_STDIO_NAMESPACE.C
// package       : Configuration (1.26)
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_STDIO_NAMESPACE.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The flag CGAL_CFG_NO_STDIO_NAMESPACE is set, if a compiler does not
//| put the IO standard library in namespace std.

#include <iosfwd>
#include <iostream>
#include <iomanip>
#include <streambuf>
#include <fstream>

// just some randomly selected symbols
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::ofstream;

int main()
{
  return 0;
}

// EOF //

