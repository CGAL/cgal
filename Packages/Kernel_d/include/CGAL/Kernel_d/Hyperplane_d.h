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
// file          : include/CGAL/Kernel_d/Hyperplane_d.h
// package       : Kernel_d
// maintainer    : Michael Seel <Michael.Seel@mpi-sb.mpg.de>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel
// coordinator   : MPI Saarbruecken (Susan.Hert@mpi-sb.mpg.de)
//
// ======================================================================
#ifndef CGAL_HYPERPLANE_D_H
#define CGAL_HYPERPLANE_D_H

CGAL_BEGIN_NAMESPACE

template <class pR>
class Hyperplane_d : public pR::Hyperplane_d_base
{ public:
  typedef typename pR::Hyperplane_d_base Base;
  typedef Hyperplane_d<pR>               Self;
  typedef pR R;
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;

  Hyperplane_d(int d=0) : Base(d) {}
  Hyperplane_d(int a, int b, int c) : Base(a,b,c) {}
  Hyperplane_d(const RT& a, const RT& b, const RT& c) : 
    Base(a,b,c) {}
  Hyperplane_d(int a, int b, int c, int d) : Base(a,b,c,d) {}
  Hyperplane_d(const RT& a, const RT& b, const RT& c, const RT& d) : 
    Base(a,b,c,d) {}

  Hyperplane_d(const Point_d<R>& p, const Direction_d<R>& dir) :
    Base(p,dir) {}

  Hyperplane_d(const Hyperplane_d<R> &h) : Base(h) {}
  Hyperplane_d(const Base& p) : Base(p) {}

#ifndef CGAL_SIMPLE_INTERFACE

  template <class InputIterator>
  Hyperplane_d(int d, InputIterator first, InputIterator last)
    : Base (d, first, last) {}

  template <class InputIterator>
  Hyperplane_d(int d, InputIterator first, InputIterator last,
               const RT& D)
    : Base (d, first, last, D) {}

  template <class ForwardIterator>
  Hyperplane_d(ForwardIterator first, ForwardIterator last, 
               const Point_d<R>& o, Oriented_side side = Oriented_side(0)) :
    Base(first,last,o,side) {}

#else
#define FIXHYPD(I) \
Hyperplane_d(int d, I first, I last) : Base (d,first,last) {} \
Hyperplane_d(int d, I first, I last, const RT& D) : Base (d,first,last,D) {}
#define FIXHYPDD(I) \
Hyperplane_d(I first, I last, const Point_d<R>& o, Oriented_side side = \
Oriented_side(0)) : Base(o.dimension()) \
{ construct_from_points(first,last,o,side); }

	     //Oriented_side(0)) : Base(first,last,o,side) {}

FIXHYPD(int*)
FIXHYPD(const int*)
FIXHYPD(RT*)
FIXHYPD(const RT*)
#undef FIXHYPD

#ifdef _MSC_VER
FIXHYPDD(typename std::vector< Point_d<R> >::iterator)
FIXHYPDD(typename std::vector< Point_d<R> >::const_iterator)
#endif // MSC

FIXHYPDD(Point_d<R>*)
FIXHYPDD(const Point_d<R>*)
#undef FIXHYPDD
#endif

  Vector_d<R> orthogonal_vector() const 
  { return Base::orthogonal_vector(); }
  Direction_d<R> orthogonal_direction() const 
  { return Base::orthogonal_direction(); }

  bool operator==(const Self& w) const
  { return Base::operator==(w); }
  bool operator!=(const Self& w) const
  { return Base::operator!=(w); }
};

CGAL_END_NAMESPACE
#endif //CGAL_HYPERPLANE_D_H
