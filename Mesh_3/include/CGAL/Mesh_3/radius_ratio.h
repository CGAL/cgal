// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

#ifndef CGAL_MESH_3_RADIUS_RATIO_H
#define CGAL_MESH_3_RADIUS_RATIO_H

namespace CGAL {

namespace Mesh_3 {
  
template <typename K>
typename K::FT
radius_ratio(const typename K::Point_3& p0,
             const typename K::Point_3& p1,
             const typename K::Point_3& p2,
             const typename K::Point_3& p3,
             K k = K())
{
  typedef typename K::FT FT;
  typename K::Compute_squared_distance_3 sq_distance =
    k.compute_squared_distance_3_object();
  typename K::Compute_squared_radius_3 comp_sq_circumradius = 
    k.compute_squared_radius_3_object();
  typename K::Compute_volume_3 volume =
    k.compute_volume_3_object();
  typename K::Compute_area_3 area = 
    k.compute_area_3_object();

  typename K::Coplanar_3 are_coplanar = 
    k.coplanar_3_object();

  if(are_coplanar(p0, p1, p2, p3))
    return FT(0.);

  const FT sq_circumradius = comp_sq_circumradius(p0,
						  p1,
						  p2,
						  p3);
  CGAL_assertion(sq_circumradius != FT(0));

  const FT triangles_area =
    area(p0, p1, p2)
    + area(p1, p2, p3)
    + area(p2, p3, p0)
    + area(p3, p0, p1);

  CGAL_assertion(triangles_area != FT(0));

  const FT cell_volume = volume(p0, p1, p2, p3);

  const FT inradius = 
    3 * CGAL::abs(cell_volume) / triangles_area;

  const FT result = FT(3) * inradius / CGAL::sqrt(sq_circumradius);
  
  CGAL_assertion(result >= FT(0));
  CGAL_assertion(result <= FT(1));

  return result;
}
  
template <typename K>
typename K::FT
radius_ratio(const typename K::Tetrahedron_3& t, K k = K())
{
  return radius_ratio(t[0], t[1], t[2], t[3], k);
}

template <typename Tetrahedron_3>
typename Kernel_traits<Tetrahedron_3>::Kernel::FT
radius_ratio(const Tetrahedron_3& t)
{  
  return radius_ratio(t, typename Kernel_traits<Tetrahedron_3>::Kernel());
}
  


} // end namespace Mesh_3
} // end namespace CGAL

#endif // end CGAL_MESH_3_RADIUS_RATIO_H
