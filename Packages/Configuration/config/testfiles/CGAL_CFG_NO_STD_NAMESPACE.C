// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_STD_NAMESPACE.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_STD_NAMESPACE.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't know the namespace std, the flag
//| CGAL_CFG_NO_STD_NAMESPACE is set. Some compilers know namespace std
//| but don't implement namespaces in general.

#include <vector.h>

using namespace std;

int main()
{
  std::vector<int> v;
  return 0;
}

// EOF //

