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
// file          : Circle_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Sven Schoenherr
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_CIRCLE_2_H
#define CGAL_CIRCLE_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Circle_2 : public R_::Circle_2_base
{
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Circle_2_base  RCircle_2;
public:
  typedef  R_   R;

    Circle_2()
      : RCircle_2() {}

    Circle_2(const CGAL::Circle_2<R> &t)
      : RCircle_2((RCircle_2&)t) {}

    Circle_2(const RCircle_2& t)
      : RCircle_2(t) {}

    Circle_2(const Point_2 &center, const FT &squared_radius,
             const Orientation &orientation)
      : RCircle_2(center, squared_radius, orientation) {}

    Circle_2(const Point_2 &center, const FT &squared_radius)
      : RCircle_2(center, squared_radius, COUNTERCLOCKWISE) {}

    Circle_2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
      : RCircle_2(p,q,r) {}

    Circle_2(const Point_2 & p, const Point_2 & q,
             const Orientation &orientation)
      : RCircle_2(p,q,orientation) {}

    Circle_2(const Point_2 & p, const Point_2 & q)
      : RCircle_2(p,q,COUNTERCLOCKWISE) {}

    Circle_2(const Point_2 & center, const Orientation& orientation)
      : RCircle_2(center,FT(0),orientation) {}

    Circle_2(const Point_2 & center)
      : RCircle_2(center,FT(0),COUNTERCLOCKWISE) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_CIRCLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Circle_2<R> &c)
{
  typedef typename R::Circle_2_base  RCircle_2;
  return os << (const RCircle_2&)c;
}

#endif // CGAL_NO_OSTREAM_INSERT_CIRCLE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_CIRCLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Circle_2<R> &c)
{
  typedef typename R::Circle_2_base  RCircle_2;
  return is >> (RCircle_2&)c;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_CIRCLE_2

CGAL_END_NAMESPACE

#endif  // CGAL_CIRCLE_2_H
