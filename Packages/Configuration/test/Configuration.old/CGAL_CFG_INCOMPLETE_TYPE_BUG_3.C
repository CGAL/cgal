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
// file          : config/testfiles/CGAL_CFG_INCOMPLETE_TYPE_BUG_3.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_INCOMPLETE_TYPE_BUG_3.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| When a class (A) refers to a not yet defined class (B), some compilers
//| give an "incomplete type error".
//| This program is used to detect a special case of this problem
//| where virtual functions and templates are involved.

#ifdef __GNUC__
#include <typeinfo>
#endif

template < class FT >
class B;

template < class FT >
class C;

template < class FT >
class A
{
public:
  virtual int f(const B<FT> &) const  = 0;
};

template < class FT >
class B : public A<FT>
{
public:
  int f(const B<FT> &) const { return 1; }
};

template < class FT >
class C : public A<FT>
{
public:
  int f(const B<FT> &) const { return 2; }
};

int main()
{
  C<int> c;

  return 0;
}

// EOF //
