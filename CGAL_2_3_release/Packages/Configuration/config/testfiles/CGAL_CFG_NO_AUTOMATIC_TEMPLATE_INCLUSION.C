// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it is able
// locate template implementation files with suffix .C  automatically.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| When template implementation files are not included in the source files,
//| a compiler may attempt to find the unincluded template bodies
//| automatically. For example, suppose that the following conditions are
//| all true.
//|
//| - template entity ABC::f is declared in file xyz.h
//| - an instantiation of ABC::f is required in a compilation
//| - no definition of ABC::f appears in the source code processed by the
//|   compilation
//| 
//| In this case, the compiler may look to see if the source file xyz.n exists,
//| where n is .c, .C, .cpp, .CPP, .cxx, .CXX, or .cc. If this feature is
//| missing, the flag CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION is set.

#include "template.h"

int main()
{
  A<float> a;
  float f = a.square(1.2);
  (void) f;

  return 0;
}

// EOF //
