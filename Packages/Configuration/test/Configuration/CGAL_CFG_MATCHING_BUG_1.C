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
// release       : $CGAL_Revision: V $
// release_date  : $CGAL_Date: 1998, July 17 $
//
// file          : config/testfiles/CGAL_CFG_MATCHING_BUG_1.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_MATCHING_BUG_1.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if the compiler does not match a function
//| argument of type typename T::A correctly.
//| (e.g. SGI 7.2)

struct Z {
  typedef int A;
};

template < class T >
void foo( T, typename T::A) 
{}

int main() {
  Z z;
  foo( z, 0);
  return 0;
}

// EOF //
