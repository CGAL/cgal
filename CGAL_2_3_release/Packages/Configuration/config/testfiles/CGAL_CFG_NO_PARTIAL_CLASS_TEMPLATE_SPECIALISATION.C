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
// file          : config/testfiles/CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION.C
// author(s)     : various
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION.C
// ---------------------------------------------------------------------

//| If a compiler doesn't support partial specialisation of class templates,
//| the flag CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION is set.

template < class T >
struct B {
  T a;
};

template < class T >
struct C {
  T b;
};

template < class T >
struct C< B < T > > {
  int c;
};

int
main()
{
  C< int > t1;
  C< B< int > > t2;
  t1.b = 1;
  t2.c = 1;
  return 0;
}

// EOF //
