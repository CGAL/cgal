// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : Triangle_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_TRIANGLE_2_H
#define CGAL_TRIANGLE_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Triangle_2 : public R_::Kernel_base::Triangle_2
{
  typedef typename R_::Point_2          Point_2;
  typedef typename R_::Kernel_base::Triangle_2  RTriangle_2;
public:
  typedef  R_                          R;

  Triangle_2()
      : RTriangle_2() {}

  Triangle_2(const CGAL::Triangle_2<R> &t)
      : RTriangle_2((RTriangle_2&)t) {}

  Triangle_2(const RTriangle_2& t)
      : RTriangle_2(t) {}

  Triangle_2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
      : RTriangle_2(p,q,r) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Triangle_2<R> &t)
{
  typedef typename  R::Kernel_base::Triangle_2  RTriangle_2;
  return os << (const RTriangle_2&)t;
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Triangle_2<R> &t)
{
  typedef typename  R::Kernel_base::Triangle_2  RTriangle_2;
  return is >> (RTriangle_2&)t;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_2

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLE_2_H
