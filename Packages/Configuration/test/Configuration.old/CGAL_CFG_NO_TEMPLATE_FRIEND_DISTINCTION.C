// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it distinguishes
// a template instantiation and a non-templated function with the same
// signature in a friend declaration.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't distinguish a template instantiation and 
//| a non-templated function with the same signature, 
//| CGAL_CFG_NO_TEMPLATE_FRIEND_DISTINCTION is set.

template < class T >
int y( const T& t) {
  return t.a;
}

template < class T >
T
operator+(const T& a, const T&) {
  return a;
}

class X {
  int a;
  friend int y<>( const X&);
  friend X operator+<>(const X&, const X&);
};

int main()
{
  X x;
  y( x);
  X z = x + x;
  return 0;
}

// EOF //

