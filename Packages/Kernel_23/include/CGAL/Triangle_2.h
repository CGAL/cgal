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
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_TRIANGLE_2_H
#define CGAL_TRIANGLE_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Aff_transformation_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Triangle_2 : public R_::Triangle_2_base
{
public:
  typedef  R_                          R;
  typedef typename R::RT               RT;
  typedef typename R::FT               FT;
  typedef typename R::Triangle_2_base  RTriangle_2;

  Triangle_2() : RTriangle_2() {}

  Triangle_2(const CGAL::Triangle_2<R> &t) : RTriangle_2((RTriangle_2&)t) {}

  Triangle_2(const RTriangle_2& t) : RTriangle_2(t) {}

  Triangle_2(const CGAL::Point_2<R> &p,
             const CGAL::Point_2<R> &q,
             const CGAL::Point_2<R> &r) : RTriangle_2(p,q,r) {}

  CGAL::Point_2<R>
  vertex(int i) const
  { return RTriangle_2::vertex(i); }

  CGAL::Point_2<R>
  operator[](int i) const
  { return vertex(i); }

  CGAL::Triangle_2<R>
  transform(const CGAL::Aff_transformation_2<R> &t) const
  { return  RTriangle_2::transform(t); }

  CGAL::Triangle_2<R>  opposite() const
  { return  CGAL::Triangle_2<R>(vertex(0), vertex(2), vertex(1)); }
};

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Triangle_2<R> &t)
{
  typedef typename  R::Triangle_2_base  RTriangle_2;
  return os << (const RTriangle_2&)t;
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Triangle_2<R> &t)
{
  typedef typename  R::Triangle_2_base  RTriangle_2;
  return is >> (RTriangle_2&)t;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_2

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLE_2_H
