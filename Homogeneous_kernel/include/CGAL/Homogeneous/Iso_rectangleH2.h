// Copyright (c) 1999
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
// Author(s)     : Stefan Schirra

#ifndef CGAL_ISO_RECTANGLEH2_H
#define CGAL_ISO_RECTANGLEH2_H

#include <CGAL/array.h>

namespace CGAL {

template <class R_>
class Iso_rectangleH2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Iso_rectangle_2      Iso_rectangle_2;

  typedef std::array<Point_2, 2>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                    R;
  typedef typename Point_2::Rep::Cartesian_coordinate_type Cartesian_coordinate_type;
  typedef typename Point_2::Rep::Homogeneous_coordinate_type Homogeneous_coordinate_type;

  Iso_rectangleH2() {}

  Iso_rectangleH2(const Point_2& p, const Point_2& q, int)
    : base(CGAL::make_array(p, q))
  {
    // I have to remove the assertions, because of Homogeneous_converter.
    // CGAL_kernel_assertion(p.x()<=q.x());
    // CGAL_kernel_assertion(p.y()<=q.y());
  }

  const Point_2 & min BOOST_PREVENT_MACRO_SUBSTITUTION () const;
  const Point_2 & max BOOST_PREVENT_MACRO_SUBSTITUTION () const;

  Bounded_side bounded_side(const Point_2& p) const;
};



template < class R >
inline
const typename Iso_rectangleH2<R>::Point_2 &
Iso_rectangleH2<R>::min BOOST_PREVENT_MACRO_SUBSTITUTION () const
{ return get_pointee_or_identity(base)[0]; }

template < class R >
inline
const typename Iso_rectangleH2<R>::Point_2 &
Iso_rectangleH2<R>::max BOOST_PREVENT_MACRO_SUBSTITUTION () const
{ return get_pointee_or_identity(base)[1]; }

template < class R >
CGAL_KERNEL_INLINE
Bounded_side
Iso_rectangleH2<R>::
bounded_side(const typename Iso_rectangleH2<R>::Point_2& p) const
{
  Oriented_side wrt_min = _where_wrt_L_wedge((this->min)(),p);
  Oriented_side wrt_max = _where_wrt_L_wedge(p,(this->max)());
  if (( wrt_min == ON_NEGATIVE_SIDE )||( wrt_max == ON_NEGATIVE_SIDE))
  {
      return ON_UNBOUNDED_SIDE;
  }
  if (  ( wrt_min == ON_ORIENTED_BOUNDARY )
      ||( wrt_max == ON_ORIENTED_BOUNDARY ) )
  {
      return ON_BOUNDARY;
  }
  return ON_BOUNDED_SIDE;
}

} //namespace CGAL

#endif // CGAL_ISO_RECTANGLEH2_H
