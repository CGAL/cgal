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
// release       : $CGAL_Revision: CGAL-0.9-I-03 $
// release_date  : $CGAL_Date: 1997/11/13 $
//
// file          : config/testfiles/CGAL_CFG_NO_STATIC_MEMBERS_IN_CLASS_TEMPLATES.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_NO_STATIC_MEMBERS_IN_CLASS_TEMPLATES.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The following flag is set if a compiler doesn't accept static members
//| in a template class.

template <class T>
class A {
  static int i;
};

template <class T>
int A<T>::i = 1;

int main() {
  A<int> a;

  return 0;
}

// EOF //
