// Copyright (c) 2008  GeometryFactory, Sophia-Antipolis (France).
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

#ifndef CGAL_MESH_3_DIHEDRAL_ANGLE_3_H
#define CGAL_MESH_3_DIHEDRAL_ANGLE_3_H

#include <CGAL/Kernel_traits.h>
#include <cmath>


namespace CGAL {
namespace Mesh_3 {
    
  
namespace details {

template <typename K>
typename K::FT
dihedral_angle_aux_compute_quotient(const typename K::Point_3& p0,
                                    const typename K::Point_3& p1,
                                    const typename K::Point_3& p2,
                                    const typename K::Point_3& p3,
                                    K k = K())
{
  typename K::Construct_triangle_3 make_triangle = 
    k.construct_triangle_3_object();
  typename K::Compute_area_3 area = 
    k.compute_area_3_object();
  typename K::Compute_squared_distance_3 sq_distance = 
    k.compute_squared_distance_3_object();
  
  return CGAL::sqrt(sq_distance(p0, p1))
      / area(make_triangle(p0, p1, p3))
      / area(make_triangle(p0, p1, p2));
}

} // end namespace details;

 
/**
 * Computes dihedral angle of planes (a,b,c) and (a,b,d)
 */
template <typename K>
typename K::FT
dihedral_angle(const typename K::Point_3& a,
               const typename K::Point_3& b,
               const typename K::Point_3& c,
               const typename K::Point_3& d,
               K k = K())
{
  typedef typename K::FT FT;
  typename K::Compute_volume_3 volume = 
    k.compute_volume_3_object();
  
  using details::dihedral_angle_aux_compute_quotient;
  
  FT quotient = dihedral_angle_aux_compute_quotient(a, b, c, d, k);
  
  return ( std::asin( FT(1.5) * volume(a, b, c, d) * quotient )
           * FT(180)
           / FT(CGAL_PI) );
}

/**
 * Computes dihedral angle of planes (a,b,c) and (a,b,d)
 */  
template <typename Point_3>
typename Kernel_traits<Point_3>::Kernel::FT
dihedral_angle(const Point_3& a, const Point_3& b,
               const Point_3& c, const Point_3& d)
{
  return dihedral_angle(a, b, c, d, typename Kernel_traits<Point_3>::Kernel());
}
  
} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_DIHEDRAL_ANGLE_3_H
