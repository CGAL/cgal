// ======================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Kernel_d/Point_d.h
// package       : Kernel_d
// maintainer    : Michael Seel <Michael.Seel@mpi-sb.mpg.de>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel
// coordinator   : MPI Saarbruecken (Susan.Hert@mpi-sb.mpg.de)
//
// ======================================================================
#ifndef CGAL_POINT_D_H
#define CGAL_POINT_D_H

CGAL_BEGIN_NAMESPACE

template <class pR>
class Point_d : public pR::Point_d_base
{ public:
  typedef typename pR::Point_d_base Base;
  typedef Point_d<pR>               Self;
  typedef pR R;
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;

  Point_d(int d=0) : Base(d) {}
  Point_d(int d, const Origin &o) : Base(d,o) {}

  Point_d(int a, int b, int c = 1) : 
    Base(RT(a),RT(b),RT(c)) {} 
  Point_d(const RT& a, const RT& b, const RT& c = 1) :
    Base(a,b,c) {}  
  Point_d(int a, int b, int c, int d) : 
    Base(RT(a),RT(b),RT(c),RT(d)) {}
  Point_d(const RT& a, const RT& b, const RT& c, const RT& d) :
    Base(a,b,c,d) {}

#ifndef CGAL_SIMPLE_INTERFACE

  template <class InputIterator>
  Point_d (int d, InputIterator first, InputIterator last)
    : Base (d, first, last) {}
  template <class InputIterator>
  Point_d(int d, InputIterator first, InputIterator last, const RT& D)
    : Base (d, first, last, D) {}

#else
#define FIXPNTD(I) \
Point_d (int d, I first, I last) : Base (d, first, last) {} \
Point_d(int d, I first, I last, const RT& D)  : Base (d, first, last, D) {}

FIXPNTD(int*)
FIXPNTD(const int*)
FIXPNTD(RT*)
FIXPNTD(const RT*)
#undef FIXPNTD
#endif

  Point_d(const Self &p) : Base(p) {}
  Point_d(const Base& p) : Base(p) {}
 
  Vector_d<R> operator-(const Origin& o) const 
  { return Base::operator-(o); }
  Vector_d<R> operator-(const Self& q) const
  { return Base::operator-(q); }
  Self operator+(const Vector_d<R>& v) const
  { return Base::operator+(v); }
  Self operator-(const Vector_d<R>& v) const
  { return Base::operator-(v); }
  Self& operator+=(const Vector_d<R>& v) 
  { return static_cast<Self&>(Base::operator+=(v)); }
  Self& operator-=(const Vector_d<R>& v)
  { return static_cast<Self&>(Base::operator-=(v)); }
  
};

CGAL_END_NAMESPACE
#endif //CGAL_POINT_D_H
