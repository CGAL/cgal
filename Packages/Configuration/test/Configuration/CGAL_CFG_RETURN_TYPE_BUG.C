// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : config/testfiles/CGAL_CFG_RETURN_TYPE_BUG.C
// package       : Configuration (2.3)
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_RETURN_TYPE_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler does not accept the overloading of a template function, when
//| the template returns a reference, while the overloading doesn't.
//| In that case, CGAL_CFG_RETURN_TYPE_BUG is set.
//| This bug shows up on VC++ 6, and VC++ 7 beta 2.

template <class T>
const T&
f( const T &t)
{
  return t;
}

double
f( const double &t)
{
  return t;
}

int main()
{
  return 0;
}
