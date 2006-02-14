// Copyright (c) 1997-2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_CIRCLE_2_H
#define CGAL_CARTESIAN_CIRCLE_2_H

#include <CGAL/utility.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_ >
class CircleC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Circle_2             Circle_2;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Triple<Point_2, FT, Orientation>         Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  CircleC2() {}

  CircleC2(const Point_2 &center, const FT &squared_radius = FT(0),
           const Orientation &orient = COUNTERCLOCKWISE) // Is this new?
  {
    CGAL_kernel_precondition( ( squared_radius >= FT(0) ) &&
                              ( orient    != COLLINEAR) );

    base = Rep(center, squared_radius, orient);
  }
 
  bool           operator==(const CircleC2 &s) const;
  bool           operator!=(const CircleC2 &s) const;

  const Point_2 & center() const
  {
    return get(base).first;
  }

  const FT & squared_radius() const
  {
    return get(base).second;
  }

  Orientation orientation() const
  {
    return get(base).third;
  }

  Circle_2 orthogonal_transform(const Aff_transformation_2 &t) const;

};




template < class R >
CGAL_KERNEL_INLINE
typename CircleC2<R>::Circle_2
CircleC2<R>::orthogonal_transform
  (const typename CircleC2<R>::Aff_transformation_2 &t) const
{
  typename R::Vector_2 vec(RT(1), RT(0) );   // unit vector // AF: was FT
  vec = vec.transform(t);             // transformed
  FT sq_scale = vec.squared_length();       // squared scaling factor

  return CircleC2<R>(t.transform(center()),
                               sq_scale * squared_radius(),
                               t.is_even() ? orientation()
                                           : CGAL::opposite(orientation()));
}



CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CIRCLE_2_H
