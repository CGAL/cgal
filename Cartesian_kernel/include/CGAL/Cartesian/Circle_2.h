// Copyright (c) 1997-2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_CIRCLE_2_H
#define CGAL_CARTESIAN_CIRCLE_2_H

#include <CGAL/Cartesian/predicates_on_points_2.h>
#include <boost/tuple/tuple.hpp>

namespace CGAL {

template <class R_ >
class CircleC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Circle_2             Circle_2;
  typedef typename R_::Point_2              Point_2;

  typedef boost::tuple<Point_2, FT, Orientation>   Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  CircleC2() {}

  explicit CircleC2(const Point_2 &center, const FT &squared_radius = FT(0),
           const Orientation &orient = COUNTERCLOCKWISE) // Is this new?
  {
    CGAL_kernel_precondition( ( squared_radius >= FT(0) ) &
                              ( orient    != COLLINEAR) );

    base = Rep(center, squared_radius, orient);
  }

  bool           operator==(const CircleC2 &s) const;
  bool           operator!=(const CircleC2 &s) const;

  const Point_2 & center() const
  {
    return get_pointee_or_identity(base).template get<0>();
  }

  const FT & squared_radius() const
  {
    return get_pointee_or_identity(base).template get<1>();
  }

  Orientation orientation() const
  {
    return get_pointee_or_identity(base).template get<2>();
  }

};

} //namespace CGAL

#endif // CGAL_CARTESIAN_CIRCLE_2_H
