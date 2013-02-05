// Copyright (c) 2000  
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
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_DIRECTION_2_H
#define CGAL_CARTESIAN_DIRECTION_2_H

#include <CGAL/array.h>
#include <CGAL/Handle_for.h>

namespace CGAL {

template < class R_ >
class DirectionC2
{
  typedef DirectionC2<R_>                   Self;
  typedef typename R_::FT                   FT;
  typedef FT                                RT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Direction_2          Direction_2;

  typedef cpp11::array<FT, 2>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:

  typedef R_                                     R;

  DirectionC2() {}
 
  DirectionC2(const FT &x, const FT &y)
    : base(CGAL::make_array(x, y)) {}

  bool operator==(const DirectionC2 &d) const;
  bool operator!=(const DirectionC2 &d) const;

  Vector_2 to_vector() const;

  const RT & dx() const
  {
      return get(base)[0];
  }
  const RT & dy() const
  {
      return get(base)[1];
  }
};

template < class R >
inline
bool
DirectionC2<R>::operator==(const DirectionC2<R> &d) const
{
  if (CGAL::identical(base, d.base))
      return true;
  return equal_direction(*this, d);
}

template < class R >
inline
bool
DirectionC2<R>::operator!=(const DirectionC2<R> &d) const
{
  return !( *this == d );
}


template < class R >
inline
typename DirectionC2<R>::Vector_2
DirectionC2<R>::to_vector() const
{
  return Vector_2(dx(), dy());
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_DIRECTION_2_H
