// Copyright (c) 2000
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

#ifndef CGAL_CARTESIAN_TRIANGLE_2_H
#define CGAL_CARTESIAN_TRIANGLE_2_H

#include <CGAL/Cartesian/predicates_on_points_2.h>
#include <CGAL/array.h>

namespace CGAL {

template <class R_>
class TriangleC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Triangle_2           Triangle_2;

  typedef std::array<Point_2, 3>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                    R;

  TriangleC2() {}

  TriangleC2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
    : base(CGAL::make_array(p, q, r)) {}


  const Point_2 &
  vertex(int i) const
  {
    if (i>2) i = i%3;
    else if (i<0) i = (i%3) + 3;
    return (i==0) ? get_pointee_or_identity(base)[0] :
      (i==1) ? get_pointee_or_identity(base)[1] :
      get_pointee_or_identity(base)[2];
  }

  const Point_2 &
  operator[](int i) const
  {
    return vertex(i);
  }

};

} //namespace CGAL

#endif // CGAL_CARTESIAN_TRIANGLE_2_H
