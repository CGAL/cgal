// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: V $
// release_date  : $CGAL_Date: 1998, July 17 $
// date          : 27 May 1998
//
// file          : config/testfiles/CGAL_CFG_INCOMPLETE_TYPE_BUG_5.C
// author(s)     : various
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_INCOMPLETE_TYPE_BUG_5.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Some compilers issue an "incomplete type error" with some STL versions
//| if a template type, which has not been instantiated yet, is put into 
//| an STL vector (e.g. SunPro 4.2 with STLport-3.01). This program is 
//| used to detect this bug.

#include <vector.h>

template <class T>
class A
{
public:
  A(const T& t) : _t(t) {}
  T _t;
};

int
main()
{
  vector<A<int> > iV;

  return 0;
}

// EOF //
