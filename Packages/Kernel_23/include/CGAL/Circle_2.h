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

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/CircleH2.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Circle_2.h>
#endif

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


  bool
  operator==(const CGAL::Circle_2<R> &t) const
  { return RCircle_2::operator==(t); }

  bool
  operator!=(const CGAL::Circle_2<R> &t) const
  { return !(*this == t); }


  CGAL::Point_2<R>
  center() const
  { return RCircle_2::center(); }

  FT
  squared_radius() const
  { return RCircle_2::squared_radius(); }

  Orientation orientation() const
  { return RCircle_2::orientation(); }


  CGAL::Circle_2<R>

  orthogonal_transform(const CGAL::Aff_transformation_2<R> &t) const
  { return  RCircle_2::orthogonal_transform(t); }

/*
  CGAL::Circle_2<R>  transform(const CGAL::Aff_transformation_2<R> &t) const
  {
    return  Circle_2::transform(t);
  }
*/

  Oriented_side
  oriented_side(const CGAL::Point_2<R> &p) const
  { return RCircle_2::oriented_side(p); }

  Bounded_side
  bounded_side(const CGAL::Point_2<R> &p) const
  { return RCircle_2::bounded_side(p); }

  bool
  has_on_boundary(const CGAL::Point_2<R> &p) const
  { return RCircle_2::has_on_boundary(p); }

  bool
  has_on_positive_side(const CGAL::Point_2<R> &p) const
  { return RCircle_2::has_on_positive_side(p); }

  bool
  has_on_negative_side(const CGAL::Point_2<R> &p) const
  { return RCircle_2::has_on_negative_side(p); }

  bool
  has_on_bounded_side(const CGAL::Point_2<R> &p) const
  { return RCircle_2::has_on_bounded_side(p); }

  bool
  has_on_unbounded_side(const CGAL::Point_2<R> &p) const
  { return RCircle_2::has_on_unbounded_side(p); }

  bool
  is_degenerate() const
  { return RCircle_2::is_degenerate(); }

  CGAL::Circle_2<R>
  opposite() const
  {
      return CGAL::Circle_2<R>(center(), squared_radius(),
                               CGAL::opposite(orientation()));
  }

  Bbox_2
  bbox() const
  { return RCircle_2::bbox(); }

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
