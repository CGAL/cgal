// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_MATCHING_BUG_2.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_MATCHING_BUG_2.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if the compiler does not match the most
//| specialized instance of a function template correctly,
//| but complains about multiple matches.
//| (e.g. SGI 7.2)

template < class T >
void foo( T, int)
{}

template < class T, class R >
void foo( T, R)
{}

int main() {
  int i;
  foo( i, 1);
  foo( i, 'a');
  return 0;
}

// EOF //
