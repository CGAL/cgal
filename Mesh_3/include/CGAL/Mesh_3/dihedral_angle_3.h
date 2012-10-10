// Copyright (c) 2008  GeometryFactory, Sophia-Antipolis (France).
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

#ifndef CGAL_MESH_3_DIHEDRAL_ANGLE_3_H
#define CGAL_MESH_3_DIHEDRAL_ANGLE_3_H

#include <CGAL/number_type_basic.h>
#include <CGAL/Kernel_traits.h>
#include <cmath>


namespace CGAL {
namespace Mesh_3 {

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
  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Construct_cross_product_vector_3 cross_product =
    k.construct_cross_product_vector_3_object();
  typename K::Compute_squared_distance_3 sq_distance =
    k.compute_squared_distance_3_object();
  typename K::Compute_scalar_product_3 scalar_product =
    k.compute_scalar_product_3_object();

  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;

  const Vector_3 ab = vector(a,b);
  const Vector_3 ac = vector(a,c);
  const Vector_3 ad = vector(a,d);

  const Vector_3 abad = cross_product(ab,ad);
  const double x = CGAL::to_double(scalar_product(cross_product(ab,ac), abad));
  const double l_ab = CGAL::sqrt(CGAL::to_double(sq_distance(a,b)));
  const double y = l_ab * CGAL::to_double(scalar_product(ac,abad));

  return FT(std::atan2(y, x) * 180 / CGAL_PI );
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
