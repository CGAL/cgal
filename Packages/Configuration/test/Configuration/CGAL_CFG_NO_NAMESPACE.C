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
// release       : $CGAL_Revision: CGAL-2.0-I-1 $
// release_date  : $CGAL_Date: 1999/02/09 $
//
// file          : config/testfiles/CGAL_CFG_NO_NAMESPACE.C
// package       : Configuration (1.21)
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_NAMESPACE.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| If a compiler doesn't know namespaces, the flag
//| CGAL_CFG_NO_NAMESPACE is set.


namespace A {
  int foo() { return 1; }
}

namespace B {
  int foo() { return 2; }
}

bool all_assertions_correct = true;

int main()
{
  all_assertions_correct &= ( A::foo() == 1);
  all_assertions_correct &= ( B::foo() == 2);
  {
    using namespace A;
    all_assertions_correct &= ( foo() == 1);
  }
  {
    using namespace B;
    all_assertions_correct &= ( foo() == 2);
  }
  return !all_assertions_correct;
}

// EOF //

