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

#include <CGAL/Point_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Circle_2 : public R_::Circle_2_base
{
public:
  typedef  R_   R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Circle_2_base  RCircle_2;

    Circle_2()
      : RCircle_2()
    {}

    Circle_2(const CGAL::Circle_2<R> &t)
      : RCircle_2((RCircle_2&)t)
    {}

    Circle_2(const RCircle_2& t)
      : RCircle_2(t)
    {}

    Circle_2(const CGAL::Point_2<R> &center,
             const FT &squared_radius,
             const Orientation &orientation)
      : RCircle_2(center, squared_radius, orientation)
    {}

    Circle_2(const CGAL::Point_2<R> &center,
             const FT &squared_radius)
      : RCircle_2(center, squared_radius, COUNTERCLOCKWISE)
    {}

    Circle_2(const CGAL::Point_2<R> &p,
             const CGAL::Point_2<R> &q,
             const CGAL::Point_2<R> &r)
      : RCircle_2(p,q,r)
    {}

    Circle_2(const CGAL::Point_2<R> & p,
             const CGAL::Point_2<R> & q,
             const Orientation &orientation)
      : RCircle_2(p,q,orientation)
    {}

    Circle_2(const CGAL::Point_2<R> & p,
             const CGAL::Point_2<R> & q)
      : RCircle_2(p,q,COUNTERCLOCKWISE)
    {}

    Circle_2(const CGAL::Point_2<R> & center,
             const Orientation& orientation)
      : RCircle_2(center,FT(0),orientation)
    {}

    Circle_2(const CGAL::Point_2<R> & center)
      : RCircle_2(center,FT(0),COUNTERCLOCKWISE)
    {}

  CGAL::Point_2<R>
  center() const
  { return RCircle_2::center(); }

  CGAL::Circle_2<R>
  orthogonal_transform(const CGAL::Aff_transformation_2<R> &t) const
  { return  RCircle_2::orthogonal_transform(t); }

  CGAL::Circle_2<R>
  opposite() const
  {
      return CGAL::Circle_2<R>(center(), squared_radius(),
                               CGAL::opposite(orientation()));
  }
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
