// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : config/testfiles/CGAL_CFG_TYPENAME_BUG.C
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_TYPENAME_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler complains about typename, when passing a dependent
//| type as template parameter, the flag CGAL_CFG_TYPENAME_BUG is set.

template < class T > struct Zap {};

struct TT { typedef int O; };

template < class T >
void foo(T) 
{
  typedef Zap< typename T::O > O;
}

int main()
{
  TT t;
  foo(t);
  return 0;
}

