// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

#ifndef CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H
#define CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 

#include <CGAL/constructions/constructions_on_weighted_points_cartesian_3.h>
#include <CGAL/predicates/predicates_on_weighted_points_cartesian_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

CGAL_BEGIN_NAMESPACE


   
//------------------ Traits class -------------------------------------

template <class R>
class Weighted_alpha_shape_euclidean_traits_3 : public 
Regular_triangulation_euclidean_traits_3<R>
{
public:
  typedef Regular_triangulation_euclidean_traits_3<R> Base;
  typedef typename Base::Compute_squared_radius_smallest_orthogonal_sphere_3
                Compute_squared_radius_smallest_orthogonal_sphere_3;
  typedef typename Base::Side_of_bounded_orthogonal_sphere_3 
                Side_of_bounded_orthogonal_sphere_3;

  //---------------------------------------------------------------------

  Compute_squared_radius_smallest_orthogonal_sphere_3 
  compute_squared_radius_3_object() const
    {
      return Compute_squared_radius_smallest_orthogonal_sphere_3();
    }
  //---------------------------------------------------------------------

  Side_of_bounded_orthogonal_sphere_3 
  side_of_bounded_sphere_3_object() const
    {
      return Side_of_bounded_orthogonal_sphere_3();
    }
};

CGAL_END_NAMESPACE

#endif //CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_3_H 
