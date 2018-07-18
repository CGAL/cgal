// Copyright (c) 2018 GeometryFactory (France).
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
// Author(s) : Maxime Gimeno
//

#ifndef CGAL_AABB_TRANSFORMED_TRAITS_H
#define CGAL_AABB_TRANSFORMED_TRAITS_H

#include <CGAL/AABB_transformed_traits_base.h>
#include <CGAL/Filtered_predicate.h>

namespace CGAL{
template<typename BaseTraits, 
         typename Kernel,
         class HasFilteredPredicates = typename Kernel::Has_filtered_predicates_tag /* = Tag_false */>
class AABB_transformed_traits
    :public AABB_transformed_traits_base<BaseTraits, Kernel> 
{};

template<typename BaseTraits, 
         typename Kernel>
class AABB_transformed_traits<BaseTraits, Kernel, Tag_true>
    :public AABB_transformed_traits_base<BaseTraits, Kernel> 
{
  typedef typename Kernel::Exact_kernel EK;
  typedef typename Kernel::Approximate_kernel FK;
  typedef typename Kernel::C2E C2E;
  typedef typename Kernel::C2F C2F;
  
  typedef Filtered_predicate<typename AABB_transformed_traits_base<BaseTraits,EK>::Split_primitives,
  typename AABB_transformed_traits_base<BaseTraits,FK>::Split_primitives,
  C2E,C2F>   Split_primitives;
  
  typedef Filtered_predicate<typename AABB_transformed_traits_base<BaseTraits,EK>::Compute_bbox,
  typename AABB_transformed_traits_base<BaseTraits,FK>::Compute_bbox,
  C2E,C2F>   Compute_bbox;
  
  typedef Filtered_predicate<typename AABB_transformed_traits_base<BaseTraits,EK>::Do_intersect,
  typename AABB_transformed_traits_base<BaseTraits,FK>::Do_intersect,
  C2E,C2F>   Do_intersect;
  
  typedef Filtered_predicate<typename AABB_transformed_traits_base<BaseTraits,EK>::Intersection,
  typename AABB_transformed_traits_base<BaseTraits,FK>::Intersection,
  C2E,C2F>   Intersection;
  
  typedef Filtered_predicate<typename AABB_transformed_traits_base<BaseTraits,EK>::Compare_distance,
  typename AABB_transformed_traits_base<BaseTraits,FK>::Compare_distance,
  C2E,C2F>   Compare_distance;
  
  typedef Filtered_predicate<typename AABB_transformed_traits_base<BaseTraits,EK>::Closest_point,
  typename AABB_transformed_traits_base<BaseTraits,FK>::Closest_point,
  C2E,C2F>   Closest_point;
  
  typedef Filtered_predicate<typename AABB_transformed_traits_base<BaseTraits,EK>::Squared_distance,
  typename AABB_transformed_traits_base<BaseTraits,FK>::Squared_distance,
  C2E,C2F>   Squared_distance;
  
  typedef Filtered_predicate<typename AABB_transformed_traits_base<BaseTraits,EK>::Equal_3,
  typename AABB_transformed_traits_base<BaseTraits,FK>::Equal_3,
  C2E,C2F>   Equal_3;
  
};

}//end CGAL
#endif // CGAL_AABB_TRANSFORMED_TRAITS_H
