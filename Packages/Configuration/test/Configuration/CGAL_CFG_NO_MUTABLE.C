// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
// date          : 21 Apr 1998
//
// file          : config/testfiles/CGAL_CFG_NO_MUTABLE.C
// author(s)     : various
// coordinator   : Utrecht University
//
// ======================================================================


// CGAL_CFG_NO_MUTABLE.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it knows
// the keyword mutable or not.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't know the keyword mutable, the flag
//| CGAL_CFG_NO_MUTABLE is set.

#include <assert.h> 

struct A {
  A() : i( 1) {}
  mutable int i;
};

int main()
{
  const A a;
  a.i = 2;
  assert( a.i == 2);
  return 0;
}

// EOF //
