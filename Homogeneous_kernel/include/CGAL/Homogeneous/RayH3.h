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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_RAYH3_H
#define CGAL_RAYH3_H

#include <utility>
#include <CGAL/kernel_config.h>
#include <CGAL/assertions.h>
#include <CGAL/Origin.h>

namespace CGAL {

template < class R_ >
class RayH3
{
   typedef typename R_::RT                   RT;
   typedef typename R_::FT                   FT;
   typedef typename R_::Point_3              Point_3;
   typedef typename R_::Line_3               Line_3;
   typedef typename R_::Direction_3          Direction_3;
   typedef typename R_::Vector_3             Vector_3;

   typedef std::pair<Point_3, Vector_3>             Rep;
   typedef typename R_::template Handle<Rep>::type  Base;

   Base base;

public:
   typedef R_                R;

    RayH3() {}

    RayH3( const Point_3& sp, const Point_3& secondp)
      : base(sp, secondp-sp) {}

    RayH3( const Point_3& sp, const Vector_3& v)
      : base(sp, v) {}

    RayH3( const Point_3& sp, const Direction_3& d)
      : base(sp, d.to_vector()) {}

    RayH3( const Point_3& sp, const Line_3& l)
      : base(sp, l.to_vector()) {}

    const Point_3 & start() const;
    const Point_3 & source() const;
    Point_3 second_point() const;
    Point_3 point(int i) const;
    Direction_3 direction() const;
    const Vector_3 & to_vector() const;
    Line_3  supporting_line() const;
    RayH3<R>   opposite() const;
    bool           has_on(const Point_3& p) const;
    bool           collinear_has_on(const Point_3 &p) const;
    bool           is_degenerate() const;

    bool           operator==(const RayH3<R>& r) const;
    bool           operator!=(const RayH3<R>& r) const;
};

template < class R >
inline
const typename RayH3<R>::Point_3 &
RayH3<R>::source() const
{ return get_pointee_or_identity(base).first; }

template < class R >
inline
const typename RayH3<R>::Point_3 &
RayH3<R>::start() const
{ return get_pointee_or_identity(base).first; }

template < class R >
inline
const typename RayH3<R>::Vector_3 &
RayH3<R>::to_vector() const
{
  return get_pointee_or_identity(base).second;
}

template < class R >
inline
typename RayH3<R>::Direction_3
RayH3<R>::direction() const
{
  return to_vector().direction();
}

template < class R >
CGAL_KERNEL_INLINE
typename RayH3<R>::Point_3
RayH3<R>::second_point() const
{ return start() + to_vector(); }

template < class R >
CGAL_KERNEL_INLINE
typename RayH3<R>::Point_3
RayH3<R>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );
  return start() + RT(i)*to_vector();
}

template < class R >
CGAL_KERNEL_INLINE
typename RayH3<R>::Line_3
RayH3<R>::supporting_line() const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return Line_3(start(), second_point() );
}

template < class R >
CGAL_KERNEL_INLINE
RayH3<R>
RayH3<R>::opposite() const
{ return RayH3<R>( start(), - direction() ); }


template < class R >
CGAL_KERNEL_INLINE
bool
RayH3<R>::has_on(const typename RayH3<R>::Point_3 &p) const
{
  return ( (  p == start() )
         ||(  Direction_3(p - start()) == direction() ) );
}

template < class R >
inline                                      /* XXX */
bool
RayH3<R>::collinear_has_on(const typename RayH3<R>::Point_3 &p) const
{ return has_on(p); }

template < class R >
inline
bool
RayH3<R>::is_degenerate() const
{ return to_vector() == NULL_VECTOR; }

template < class R >
CGAL_KERNEL_INLINE
bool
RayH3<R>::operator==(const RayH3<R>& r) const
{ return ( (start() == r.start() )&&( direction() == r.direction() ) ); }

template < class R >
CGAL_KERNEL_INLINE
bool
RayH3<R>::operator!=( const RayH3<R>& r) const
{ return !operator==(r); }

} //namespace CGAL

#endif // CGAL_RAYH3_H
