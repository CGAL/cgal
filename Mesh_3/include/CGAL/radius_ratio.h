// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_RADIUS_RATIO_H
#define CGAL_RADIUS_RATIO_H

#include <CGAL/Point_3.h>

namespace CGAL {

template <typename K>
typename K::FT
radius_ratio(const Point_3<K>& p0,
             const Point_3<K>& p1,
             const Point_3<K>& p2,
             const Point_3<K>& p3)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;
  typename K::Compute_squared_distance_3 sq_distance =
    K().compute_squared_distance_3_object();
  typename K::Construct_circumcenter_3 circumcenter = 
    K().construct_circumcenter_3_object();
  typename K::Compute_volume_3 volume =
    K().compute_volume_3_object();
  typename K::Compute_area_3 area = 
    K().compute_area_3_object();

  const Point_3 center = circumcenter(p0,
                                      p1,
                                      p2,
                                      p3);

  const FT sq_circumradius = sq_distance(center, p0);

  CGAL_assertion(sq_circumradius != FT(0));

  FT triangles_area = FT(0);
  triangles_area += area(p0, p1, p2);
  triangles_area += area(p1, p2, p3);
  triangles_area += area(p2, p3, p0);
  triangles_area += area(p3, p0, p1);

  CGAL_assertion(triangles_area != FT(0));

  const FT cell_volume = volume(p0,
                                p1,
                                p2,
                                p3);

  const FT inradius = 
    3 * CGAL::abs(cell_volume) / triangles_area;

  const FT result = FT(3) * inradius / CGAL::sqrt(sq_circumradius);
  
  CGAL_assertion(result >= FT(0));
  CGAL_assertion(result <= FT(1));

  return result;
}

} // end namespace CGAL

#endif // end CGAL_RADIUS_RATIO_H
