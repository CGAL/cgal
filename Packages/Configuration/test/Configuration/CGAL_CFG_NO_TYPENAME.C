// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_TYPENAME.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================


// CGAL_CFG_NO_TYPENAME.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it knows
// the keyword typename or not.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't know the keyword typename, the flag
//| CGAL_CFG_NO_TYPENAME is set.

#include <cassert>

struct X {
  typedef int A;
};

template < class T >
struct Y {
  typedef typename T::A A;
  Y( const A& a) : i( a) {}
  A i;
};

int main()
{
  Y< X > i( 1);
  assert( i.i == 1);
  return 0;
}

// EOF //

