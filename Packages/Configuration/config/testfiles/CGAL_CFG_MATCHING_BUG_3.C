// ======================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.4-I-8 $
// release_date  : $CGAL_Date: 2001/09/11 $
//
// file          : config/testfiles/CGAL_CFG_MATCHING_BUG_3.C
// package       : Configuration (2.12)
// maintainer    : Geert-Jan Giezeman <geert@cs.uu.nl>
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_MATCHING_BUG_3.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if the compiler does not match function arguments 
//| of pointer type correctly, when the return type depends on 
//| the parameter's type. (e.g. sun C++ 5.3)

template < class T > struct A     { typedef typename T::CCC CCC; };
template < class T > struct A<T*> { typedef typename T::CCC CCC; };

template < class T >
typename A< T >::CCC
foo(T)
{
  typedef typename A< T >::CCC C;
  return C();
}

struct B { typedef int CCC; };

int main()
{
  B *p = 0;
  return foo(p);
}

// EOF //
