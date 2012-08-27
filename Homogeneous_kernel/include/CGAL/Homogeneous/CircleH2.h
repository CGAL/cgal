// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Sven Schoenherr
//                 Stefan Schirra

#ifndef CGAL_CIRCLEH2_H
#define CGAL_CIRCLEH2_H

#include <CGAL/Interval_nt.h>
#include <boost/tuple/tuple.hpp>

namespace CGAL {

template <class R_>
class CircleH2
{
    typedef typename R_::FT                   FT;
    typedef typename R_::RT                   RT;
    typedef typename R_::Point_2              Point_2;

    typedef boost::tuple<Point_2, FT, Orientation>   Rep;
    typedef typename R_::template Handle<Rep>::type  Base;

    Base base;

public:
    typedef R_                                    R;

    CircleH2() {}

    CircleH2(const Point_2& p, const Point_2& q, const Point_2& r)
    {
      Orientation o = CGAL::orientation( p, q, r);
      CGAL_kernel_precondition( o != COLLINEAR);

      Point_2    cp   = circumcenter( p, q, r);
      FT         sq_r = squared_distance( p, cp);

      base = Rep(cp, sq_r, o);
    }

    CircleH2(const Point_2& p, const Point_2& q, const Orientation& o)
    {
      CGAL_kernel_precondition( o != COLLINEAR);

      if ( p != q)
      {
         Point_2    cp   = midpoint( p, q);
         FT         sq_r = squared_distance( cp, p);
         base = Rep(cp, sq_r, o);
      }
      else
         base = Rep(p, FT( 0), o);
    }

    CircleH2(const Point_2& cp, const FT& squared_radius,
             const Orientation& o)
    {
      CGAL_precondition( ( ! CGAL_NTS is_negative( squared_radius)) &&
                         ( o != COLLINEAR ) );
      base = Rep(cp, squared_radius, o);
    }

    const Point_2 &
    center() const;

    Orientation
    orientation() const;

    const FT &
    squared_radius() const;

    CircleH2<R>
    opposite() const;

    Oriented_side
    oriented_side(const Point_2& ) const;

    Bounded_side
    bounded_side(const Point_2& ) const;

    bool  operator==( const CircleH2<R>& ) const;
    bool  operator!=( const CircleH2<R>& ) const;
    bool  has_on_positive_side(const Point_2& ) const;
    bool  has_on_negative_side(const Point_2& ) const;
    bool  has_on_boundary( const Point_2& ) const;
    bool  has_on_bounded_side( const Point_2& ) const;
    bool  has_on_unbounded_side(const Point_2&) const;
    bool  is_degenerate() const;

    // bool  oriented_equal( const CircleH2<R>& ) const;
    // bool  unoriented_equal( const CircleH2<R>& ) const;
};

template <class R>
inline
const typename CircleH2<R>::Point_2 &
CircleH2<R>::center() const
{ return get(base).template get<0>(); }

template <class R>
inline
const typename CircleH2<R>::FT &
CircleH2<R>::squared_radius() const
{ return get(base).template get<1>(); }

template <class R>
CGAL_KERNEL_INLINE
CircleH2<R>
CircleH2<R>::opposite() const
{
  return CircleH2<R>( center(),
                          squared_radius(),
                          CGAL::opposite( orientation() ) );
}

template <class R>
inline
Orientation
CircleH2<R>::orientation() const
{ return get(base).template get<2>(); }

template <class R>
CGAL_KERNEL_INLINE
Oriented_side
CircleH2<R>::oriented_side( const typename CircleH2<R>::Point_2& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  Comparison_result vgl = CGAL_NTS compare( sq_dist, sq_rad );
  Oriented_side rel_pos = (vgl == LARGER ) ?
                                   ON_NEGATIVE_SIDE :
                                   ( (vgl == SMALLER ) ?
                                          ON_POSITIVE_SIDE :
                                          ON_ORIENTED_BOUNDARY);
  if (orientation() == POSITIVE)
  { return rel_pos; }
  else       // NEGATIVE
  { return CGAL::opposite( rel_pos ); }
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_positive_side(const typename CircleH2<R>::Point_2& p) const
{
  if ( orientation() == POSITIVE )
  { return (has_on_bounded_side(p) ); }
  else
  { return (has_on_unbounded_side(p) ); }
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_boundary(const typename CircleH2<R>::Point_2& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_dist == sq_rad );
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_negative_side( const typename CircleH2<R>::Point_2&p) const
{
  if ( orientation() == NEGATIVE )
  {
      return (has_on_bounded_side(p) );
  }
  else
  {
      return (has_on_unbounded_side(p) );
  }
}

template <class R>
CGAL_KERNEL_INLINE
Bounded_side
CircleH2<R>::bounded_side(const typename CircleH2<R>::Point_2& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  Comparison_result vgl = CGAL_NTS compare( sq_dist, sq_rad );
  return  (vgl == LARGER ) ? ON_UNBOUNDED_SIDE :
                                   ( (vgl == SMALLER ) ?
                                          ON_BOUNDED_SIDE :
                                          ON_BOUNDARY);
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_bounded_side(const typename CircleH2<R>::Point_2& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_dist < sq_rad );
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_unbounded_side(const typename CircleH2<R>::Point_2&p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_rad < sq_dist );
}

template <class R>
inline
bool
CircleH2<R>::is_degenerate() const
{ return ( squared_radius() == FT(0) ); }

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::operator==(const CircleH2<R>& c) const
{
  return  ( center() == c.center() )
        &&( squared_radius() == c.squared_radius() )
        &&( orientation() == c.orientation() );
}

template <class R>
inline
bool
CircleH2<R>::operator!=(const CircleH2<R>& c) const
{ return !(*this == c); }

} //namespace CGAL

#endif // CGAL_CIRCLEH2_H
