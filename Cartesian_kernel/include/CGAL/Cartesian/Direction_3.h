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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_DIRECTION_3_H
#define CGAL_CARTESIAN_DIRECTION_3_H

#include <CGAL/array.h>
#include <CGAL/Handle_for.h>
#include <CGAL/predicates/kernel_ftC3.h>

namespace CGAL {

template < class R_ >
class DirectionC3
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Direction_3          Direction_3;

  typedef cpp11::array<FT, 3>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:

  typedef R_                                R;

  DirectionC3() {}

  explicit DirectionC3(const Vector_3 &v)
    : base(CGAL::make_array(v.x(), v.y(), v.z())) {}
  // { *this = v.direction(); }

  explicit DirectionC3(const Line_3 &l)
  { *this = l.rep().direction(); }

  explicit DirectionC3(const Ray_3 &r)
  { *this = r.direction(); }

  explicit DirectionC3(const Segment_3 &s)
  { *this = s.direction(); }

  DirectionC3(const FT &x, const FT &y, const FT &z)
    : base(CGAL::make_array(x, y, z)) {}

  typename R::Boolean   operator==(const DirectionC3 &d) const;
  typename R::Boolean   operator!=(const DirectionC3 &d) const;

  Vector_3       to_vector() const;
  Vector_3       vector() const { return to_vector(); }

  const FT & dx() const
  {
      return get_pointee_or_identity(base)[0];
  }
  const FT & dy() const
  {
      return get_pointee_or_identity(base)[1];
  }
  const FT & dz() const
  {
      return get_pointee_or_identity(base)[2];
  }

  const FT & hdx() const
  {
      return dx();
  }
  const FT & hdy() const
  {
      return dy();
  }
  const FT & hdz() const
  {
      return dz();
  }
  FT hw() const
  {
      return FT(1);
  }
};

template < class R >
inline
typename R::Boolean
DirectionC3<R>::operator==(const DirectionC3<R> &d) const
{
  if (CGAL::identical(base, d.base))
      return true;
  return equal_directionC3(dx(), dy(), dz(), d.dx(), d.dy(), d.dz());
}

template < class R >
inline
typename R::Boolean
DirectionC3<R>::operator!=(const DirectionC3<R> &d) const
{
  return !(*this == d);
}

template < class R >
inline
typename DirectionC3<R>::Vector_3
DirectionC3<R>::to_vector() const
{
  return Vector_3(dx(), dy(), dz());
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_DIRECTION_3_H
