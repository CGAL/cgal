// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
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
// file          : config/testfiles/CGAL_CFG_VC7_PRIVATE_TYPE_BUG.C
// package       : Configuration (2.3)
// author(s)     : Sylvain Pion
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_VC7_PRIVATE_TYPE_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This is a test-case for a bug in VC++ 7.0 beta2 that occurs in the kernel.
//| When the bug is present, CGAL_CFG_VC7_PRIVATE_TYPE_BUG is set.

template <class R>
class A
{
  typedef typename R::B B;
public:
 
  A() :i(0) {}
 
  int i;
};
 
template <class R>
class B
: public R::A
{
public:
  typedef typename R::A base;
  B() : base() {}
};
 
template <class FT>
class RR
{
public:
  typedef ::A<RR> A;
  typedef ::B<RR> B;
};
 
int main()
{
  typedef RR<int> r;
  A<r> a;
  B<r> b;
  (void) a;
  (void) b;
  return 0;
}
