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
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_SEGMENT_2_H
#define CGAL_CARTESIAN_SEGMENT_2_H

#include <CGAL/array.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

namespace CGAL {
template < class R_ >
class SegmentC2
{
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Segment_2            Segment_2;

  typedef cpp11::array<Point_2, 2>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  SegmentC2() 
  {}

  SegmentC2(const Point_2 &sp, const Point_2 &ep)
    : base(CGAL::make_array(sp, ep))
  {}

  const Point_2 &   
  source() const
  {
      return get_pointee_or_identity(base)[0];
  }
  
  const Point_2 &   
  target() const
  {
      return get_pointee_or_identity(base)[1];
  }

};

} //namespace CGAL

#endif // CGAL_CARTESIAN_SEGMENT_2_H
