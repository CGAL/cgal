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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb

#ifndef CGAL_MESH_3_MIN_DIHEDRAL_ANGLE_H
#define CGAL_MESH_3_MIN_DIHEDRAL_ANGLE_H

#include <CGAL/license/Mesh_3.h>


#include <cmath>
#include <CGAL/Kernel_traits.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

namespace CGAL {

namespace Mesh_3 {

#ifdef CGAL_MESH_3_OLD_MINIMUM_DIHEDRAL_ANGLE 

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

  typename K::Compute_squared_distance_3 sq_distance = 
    k.compute_squared_distance_3_object();
  typename K::Compute_volume_3 volume = 
    k.compute_volume_3_object();
  typename K::Compute_area_3 area = 
    k.compute_area_3_object();

  FT a_012 = area(p0,p1,p2);
  FT a_013 = area(p0,p1,p3);
  FT a_123 = area(p1,p2,p3);
  FT a_023 = area(p0,p2,p3);

  FT min_quotient = 
    CGAL::sqrt(sq_distance(p0, p1)) / a_012 / a_013;

  min_quotient = 
    (CGAL::min)(min_quotient,
                CGAL::sqrt(sq_distance(p0, p2)) / a_012 / a_023);
  min_quotient = 
    (CGAL::min)(min_quotient,
               CGAL::sqrt(sq_distance(p0, p3)) / a_013 / a_023);
  min_quotient = 
    (CGAL::min)(min_quotient,
               CGAL::sqrt(sq_distance(p1, p2)) / a_012 / a_123);
  min_quotient = 
    (CGAL::min)(min_quotient,
               CGAL::sqrt(sq_distance(p1, p3)) / a_013 / a_123);
  min_quotient = 
    (CGAL::min)(min_quotient,
               CGAL::sqrt(sq_distance(p2, p3)) / a_023 / a_123);

  const FT result (std::asin( FT(1.5) * volume(p0, p1, p2, p3) * min_quotient )
    * FT(180) / FT(CGAL_PI));
  
  return CGAL::abs(result);
}
  
#else // not CGAL_MESH_3_OLD_MINIMUM_DIHEDRAL_ANGLE 

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

  typename K::Compute_determinant_3 determinant = 
    k.compute_determinant_3_object();
  typename K::Construct_cross_product_vector_3 cp =
    k.construct_cross_product_vector_3_object();
  typename K::Compute_scalar_product_3 sp =
    k.compute_scalar_product_3_object();
  typename K::Construct_vector_3 cv =
    k.construct_vector_3_object();

  typename K::Vector_3 v01 = cv(p0,p1);
  typename K::Vector_3 v02 = cv(p0,p2);
  typename K::Vector_3 v03 = cv(p0,p3);
  typename K::Vector_3 v12 = cv(p1,p2);
  typename K::Vector_3 v13 = cv(p1,p3);
  typename K::Vector_3 v23 = cv(p2,p3);

  typename K::Vector_3 v_01_02 = cp(v01,v02);
  FT a_012 = v_01_02*v_01_02;


  typename K::Vector_3 v_01_03 = cp(v01,v03);
  FT a_013 = v_01_03*v_01_03;


  typename K::Vector_3 v_12_13 = cp(v12,v13);
  FT a_123 = v_12_13*v_12_13;

  typename K::Vector_3 v_02_03 = cp(v02,v03);
  FT a_023 = v_02_03*v_02_03;

  FT min_quotient = sp(v01,v01) / (a_012 * a_013);
  min_quotient = (CGAL::min)(min_quotient,
                             sp(v02,v02) / (a_012 * a_023));
  min_quotient = (CGAL::min)(min_quotient,
                             (v03*v03) / (a_013 * a_023));
  min_quotient = (CGAL::min)(min_quotient,
                             sp(v12,v12) / (a_012 * a_123));
  min_quotient = (CGAL::min)(min_quotient,
                             sp(v13,v13) / (a_013 * a_123));
  min_quotient = (CGAL::min)(min_quotient,
                             sp(v23,v23) / (a_023 * a_123));
  min_quotient =  sqrt(min_quotient);

  return CGAL::abs(std::asin(determinant(v01, v02, v03) * min_quotient)
                   * FT(180) / FT(CGAL_PI));
}

#endif //  CGAL_MESH_3_OLD_MINIMUM_DIHEDRAL_ANGLE 

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

template<typename C3T3>
typename C3T3::Triangulation::Geom_traits::FT
minimum_dihedral_angle_in_c3t3(const C3T3& c3t3)
{
  typedef typename C3T3::Triangulation::Geom_traits K;
  typename K::FT min_angle = (typename K::FT)(90.);

  typename C3T3::Cells_in_complex_iterator cit;
  for(cit = c3t3.cells_in_complex_begin();
      cit != c3t3.cells_in_complex_end();
      ++cit)
  {
    min_angle = (std::min)(min_angle, 
      minimum_dihedral_angle(c3t3.triangulation().tetrahedron(cit)));
  }
  return min_angle;
}
  

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_MIN_DIHEDRAL_ANGLE_H
