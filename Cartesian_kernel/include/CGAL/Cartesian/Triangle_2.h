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

  typedef cpp11::array<Point_2, 3>          Rep;
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
    return (i==0) ? get(base)[0] :
      (i==1) ? get(base)[1] :
      get(base)[2];
  }
  
  const Point_2 &
  operator[](int i) const
  {
    return vertex(i);
  }

};

} //namespace CGAL

#endif // CGAL_CARTESIAN_TRIANGLE_2_H
