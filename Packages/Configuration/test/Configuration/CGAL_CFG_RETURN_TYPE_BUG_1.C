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
// release_date  : $CGAL_Date: 1997/12/01 $
//
// file          : config/testfiles/CGAL_CFG_RETURN_TYPE_BUG_1.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_RETURN_TYPE_BUG_1.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The following flag is set if a compiler does not allow the return
//| type T::A of a function or a class method.

template <class T>
typename T::A f(T) { return T::A(); }

template <class T> class B;

template <class T>
typename T::A g(B<T>) { return T::A(); }

int main()
{
  return 0;
}

// EOF //
