// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three kinds:
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the CGAL Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the CGAL Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Max-Planck-Institute Saarbrucken
// (Germany), RISC Linz (Austria), and Tel-Aviv University (Israel).
//
// ============================================================================
//
// release       : CGAL-1.0a
// date          : 27 May 1998
//
// file          : config/testfiles/CGAL_CFG_INCOMPLETE_TYPE_BUG_5.C
// author(s)     : various
//
// email         : cgal@cs.uu.nl
//
// ============================================================================

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
