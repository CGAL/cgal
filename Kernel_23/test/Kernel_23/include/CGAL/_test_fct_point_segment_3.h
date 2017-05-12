// Copyright (c) 2001  
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
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri
 

#ifndef CGAL__TEST_FCT_POINT_3_SEGMENT_3_H
#define CGAL__TEST_FCT_POINT_3_SEGMENT_3_H


template <class R>
bool
_test_fct_point_segment_3(const R& )
{
  std::cout << "Testing functions Point_3 Segment_3" ;
  R k;
  typename R::Construct_projected_point_and_location_3 project = k.construct_projected_point_and_location_3_object();
  CGAL::Triangle_3<R> t(CGAL::Point_3<R>(0, 0, 0), CGAL::Point_3<R>( 3, 0, 0), CGAL::Point_3<R>( 0, 3, 0));

  typename R::Construct_projected_point_and_location_3::result_type pdi;

  pdi = project(t, CGAL::Point_3<R>(-1, -1, 0));
  assert((pdi.dimension == 0) && (pdi.index == 0));

  pdi = project(t, CGAL::Point_3<R>(4, -1, 0));
  assert((pdi.dimension == 0) && (pdi.index == 1));

  pdi = project(t, CGAL::Point_3<R>(-1, 4, 0));
  assert((pdi.dimension == 0) && (pdi.index == 2));

  pdi = project(t, CGAL::Point_3<R>(-1, 1, 0));
  assert((pdi.dimension == 1) && (pdi.index == 1));
  
  pdi = project(t, CGAL::Point_3<R>(1, -1, 0));
  assert((pdi.dimension == 1) && (pdi.index == 2));
  
  pdi = project(t, CGAL::Point_3<R>(3, 3, 0));
  assert((pdi.dimension == 1) && (pdi.index == 0));
  
  pdi = project(t, CGAL::Point_3<R>(1, 1, 1));
  assert(pdi.dimension == 2);

 std::cout << "done" << std::endl;
 return true;
}
#endif // CGAL__TEST_FCT_POINT_3_SEGMENT_3_H
