// Copyright (c) 1999,2004   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion
//                 Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H

#include <CGAL/basic.h>

namespace CGAL {



template < class K, class Weight = typename K::RT >
class Regular_triangulation_euclidean_traits_3
  : public K
{
public:
  typedef K                                          Kernel;
  typedef typename K::FT                             FT;
  typedef typename K::Point_3                        Bare_point;
  typedef typename K::Weighted_point_3               Weighted_point;
  typedef Weighted_point                             Weighted_point_3;
  typedef Weighted_point                             Point_3;

  typedef Regular_triangulation_euclidean_traits_3<K, Weight> Self;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3 Point;

  typedef typename K::Power_side_of_oriented_power_sphere_3                 Power_side_of_oriented_power_sphere_3;
  typedef typename K::Compare_power_distance_3     Compare_power_distance_3;
  typedef typename K::Construct_weighted_circumcenter_3 Construct_weighted_circumcenter_3;


  typedef typename K::Power_side_of_bounded_power_sphere_3
                                Power_side_of_bounded_power_sphere_3;
  typedef typename K::Side_of_bounded_orthogonal_sphere_3 
                                Side_of_bounded_orthogonal_sphere_3;
  typedef typename K::Compute_squared_radius_smallest_orthogonal_sphere_3 
                Compute_squared_radius_smallest_orthogonal_sphere_3;
  typedef typename K::Compute_power_product_3     Compute_power_product_3;
  typedef typename K::Compute_power_distance_to_power_sphere_3 
                                       Compute_power_distance_to_power_sphere_3;
  typedef typename K::Compare_weighted_squared_radius_3 
                                       Compare_weighted_squared_radius_3;


  Power_side_of_oriented_power_sphere_3   power_side_of_oriented_power_sphere_3_object() const
  { return K().power_side_of_oriented_power_sphere_3_object(); }

  Compare_power_distance_3 compare_power_distance_3_object() const
  { return K().compare_power_distance_3_object(); }

  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object() const
  { return K().construct_weighted_circumcenter_3_object(); }


  Power_side_of_bounded_power_sphere_3
  power_side_of_bounded_power_sphere_3_object() const
  { return K().in_smallest_orthogonal_sphere_3_object(); }

  Side_of_bounded_orthogonal_sphere_3
  side_of_bounded_orthogonal_sphere_3_object() const
  { return K().side_of_bounded_orthogonal_sphere_3_object(); }

  Compute_power_product_3
  compute_power_product_3_object() const
  { return K().compute_power_product_3_object(); }

  Compute_squared_radius_smallest_orthogonal_sphere_3
  compute_squared_radius_smallest_orthogonal_sphere_3_object() const
  { return K().compute_squared_radius_smallest_orthogonal_sphere_3_object(); }

  Compute_power_distance_to_power_sphere_3
  compute_power_distance_to_power_sphere_3_object() const
  {return  K().compute_power_distance_to_power_sphere_3_object(); }
  
  Compare_weighted_squared_radius_3
  compare_weighted_squared_radius_3_object() const
  {return K().compare_weighted_squared_radius_3_object();  }

};


} //namespace CGAL



#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_3_H
