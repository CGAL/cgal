// Copyright (c) 2007-2009  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU, Stephane Tayeb

#ifndef CGAL_MESH_3_MIN_DIHEDRAL_ANGLE_H
#define CGAL_MESH_3_MIN_DIHEDRAL_ANGLE_H

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <cmath>

namespace CGAL {

namespace Mesh_3 {
  
  namespace details {
    
    template <typename K>
    typename K::FT
    min_dihedral_angle_aux_compute_quotient(const typename K::Point_3& p0,
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
  
  
template <typename K>
typename K::FT
minimum_dihedral_angle(
     const typename K::Point_3& p0,
     const typename K::Point_3& p1,
     const typename K::Point_3& p2,
     const typename K::Point_3& p3,
     K k = K())
{
  typedef typename K::FT FT;
  typename K::Compute_volume_3 volume = 
    k.compute_volume_3_object();

  using details::min_dihedral_angle_aux_compute_quotient;

  FT min_quotient = 
    min_dihedral_angle_aux_compute_quotient(p0, p1, p2, p3, k);

  min_quotient = 
    (std::min)(min_quotient,
       min_dihedral_angle_aux_compute_quotient(p0, p2, p1, p3, k));
  min_quotient = 
    (std::min)(min_quotient,
       min_dihedral_angle_aux_compute_quotient(p0, p3, p1, p2, k));
  min_quotient = 
    (std::min)(min_quotient,
       min_dihedral_angle_aux_compute_quotient(p1, p2, p0, p3, k));
  min_quotient = 
    (std::min)(min_quotient,
       min_dihedral_angle_aux_compute_quotient(p1, p3, p0, p2, k));
  min_quotient = 
    (std::min)(min_quotient,
       min_dihedral_angle_aux_compute_quotient(p2, p3, p0, p1, k));

  const FT result (std::asin( FT(1.5) * volume(p0, p1, p2, p3) * min_quotient )
    * FT(180) / FT(CGAL_PI));
  
  return CGAL::abs(result);
}
  
  
template <typename K>
typename K::FT
minimum_dihedral_angle(const typename K::Tetrahedron_3& t, K k = K() )
{
  return minimum_dihedral_angle(t[0],t[1],t[2],t[3],k);
}

template <typename Tetrahedron_3>
typename Kernel_traits<Tetrahedron_3>::Kernel::FT
minimum_dihedral_angle(const Tetrahedron_3& t )
{
  return minimum_dihedral_angle(t, typename Kernel_traits<Tetrahedron_3>::Kernel() );
}
  

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_MIN_DIHEDRAL_ANGLE_H
