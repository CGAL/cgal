// Copyright (c) 1999,2016
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
// Author(s)     : Mariette Yvinec
//                 Sylvain Pion

#ifndef CGAL_HOMOGENEOUS_WEIGHTED_POINT_3_H
#define CGAL_HOMOGENEOUS_WEIGHTED_POINT_3_H

#include <CGAL/Handle_for.h>
#include <CGAL/Origin.h>

#include <boost/tuple/tuple.hpp>

#include <iostream>

namespace CGAL {

template < class R_ >
class Weighted_pointH3
{
  typedef typename R_::Point_3                     Point_3;
  typedef typename R_::FT                          FT;
  typedef FT                                       Weight;

  typedef boost::tuple<Point_3, Weight>            Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  Weighted_pointH3 ()
  {}

  Weighted_pointH3(const Origin &o)
    : base(o,0) {}

  explicit
  Weighted_pointH3 (const Point_3 &p)
    : base(p,0)
  {}

  Weighted_pointH3 (const Point_3 &p, const Weight &w)
    : base(p,w)
  {}

  // Constructors from coordinates are also provided for convenience, except
  // that they are only from Cartesian coordinates, and with no weight, so as
  // to avoid any potential ambiguity between the homogeneous weight and the
  // power weight (it should be easy enough to pass a Point_3 explicitly in those
  // cases).

  Weighted_pointH3 (const FT &x, const FT &y, const FT &z)
    : base(Point_3(x, y, z), 0)
  {}

  const Point_3 & point() const
  {
    return get_pointee_or_identity(base).template get<0>();
  }

  const Weight & weight() const
  {
    return get_pointee_or_identity(base).template get<1>();
  }
};

template < class R_ >
std::ostream &
operator<<(std::ostream &os, const Weighted_pointH3<R_> &p)
{
  switch(get_mode(os))
  {
  case IO::ASCII :
    return os << p.point() <<  " " << p.weight();
  case IO::BINARY :
    os << p.point();
    write(os, p.weight());
    return os;
  default:
    return os << "Weighted_point_3(" << p.point() << ", " << p.weight() << ")";
  }
}

template < class R_ >
std::istream &
operator>>(std::istream &is, Weighted_pointH3<R_> &wp)
{
  typename Weighted_pointH3<R_>::Weight w;
  typename Weighted_pointH3<R_>::Point_3 p;
  is >> p;
  if(!is) return is;
  if(is_ascii(is))
    is >> w;
  else
    read(is, w);
  if (is)
    wp = Weighted_point_3<R_>(p,w);
  return is;
}

} // namespace CGAL

#endif // CGAL_HOMOGENEOUS_WEIGHTED_POINT_3_H
