
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

#ifndef CGAL_AABB_TRANSFORMED_TRAITS_BASE_H
#define CGAL_AABB_TRANSFORMED_TRAITS_BASE_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Default.h>
#include <CGAL/intersections.h>
#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <CGAL/internal/AABB_tree/Is_ray_intersection_geomtraits.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>

#include <CGAL/Aff_transformation_3.h>
#include <boost/optional.hpp>
#include <boost/bind.hpp>

/// \file AABB_transformed_traits_base.h

namespace CGAL {
// forward declaration
template< typename AABBTraits>
class AABB_tree;
/// \addtogroup PkgAABB_tree
/// @{

/// \tparam BaseTraits a model of `CGAL::AABBTraits`
/// 
/// \sa `AABBTraits`
/// \sa `AABB_tree`
/// \sa `AABBPrimitive`
/// \sa `AABBPrimitiveWithSharedData`

  template<typename BaseTraits, 
           typename Kernel>
class AABB_transformed_traits_base:
  public BaseTraits
{
public:

  //Constructor
  AABB_transformed_traits_base(const Aff_transformation_3<Kernel>& transf = Aff_transformation_3<Kernel>(IDENTITY))
    :m_transfo(transf)
  {}
  // AABBTraits concept types
  typedef typename BaseTraits::FT FT;
  typedef typename BaseTraits::Point_3 Point_3;
  typedef typename BaseTraits::Primitive Primitive;
  typedef typename BaseTraits::Bounding_box Bounding_box;
  typedef typename BaseTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename BaseTraits::Object_and_primitive_id  Object_and_primitive_id;
  template<typename Query>
  struct Intersection_and_primitive_id {
    typedef typename BaseTraits::template Intersection_and_primitive_id<Query>::Intersection_type Intersection_type;

    typedef typename BaseTraits::template Intersection_and_primitive_id<Query>::Type Type;
  };
  
  //SearchGeomTriats_3 concept types
  typedef typename BaseTraits::Iso_cuboid_3  Iso_cuboid_3;
  typedef typename BaseTraits::Sphere_3  Sphere_3;
  typedef typename BaseTraits::Construct_iso_cuboid_3  Construct_iso_cuboid_3;
  typedef typename BaseTraits::Construct_min_vertex_3  Construct_min_vertex_3;
  typedef typename BaseTraits::Construct_max_vertex_3  Construct_max_vertex_3;
  typedef typename BaseTraits::Construct_center_3  Construct_center_3;
  typedef typename BaseTraits::Compute_squared_radius_3  Compute_squared_radius_3;
  typedef typename BaseTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_3;
  typedef typename BaseTraits::Construct_cartesian_const_iterator_3  Construct_cartesian_const_iterator_3;
  
  //Splitting
    typedef typename BaseTraits::Split_primitives            Split_primitives;
  
  typedef typename  BaseTraits::Compute_bbox                 Compute_bbox;
  
  //Intersections
  class Do_intersect {
    const AABB_transformed_traits_base<BaseTraits, Kernel>& m_traits;
  public:
    Do_intersect(const AABB_transformed_traits_base<BaseTraits, Kernel>& traits)
      :m_traits(traits) {}

    template<typename Query>
    bool operator()(const Query& q, const Bounding_box& bbox) const
    {
      Point_3 min(bbox.xmin(), bbox.ymin(), bbox.zmin()),
      max(bbox.xmax(), bbox.ymax(), bbox.zmax());
      
      min = m_traits.transformation().transform(min);
      max = m_traits.transformation().transform(max);
      Bounding_box transfo_box(to_double(min.x()), to_double(min.y()), to_double(min.z()),
                               to_double(max.x()), to_double(max.y()), to_double(max.z()));
      bool res = CGAL::do_intersect(q, transfo_box);
      return res;
    }

    template<typename Query>
    bool operator()(const Query& q, const Primitive& pr) const
    {
      return Kernel().do_intersect_3_object()(q, internal::Primitive_helper<BaseTraits>::get_datum(pr,m_traits));
    }
    
    // intersection with AABB-tree
    template<typename AABBTraits>
    bool operator()(const CGAL::AABB_tree<AABBTraits>& other_tree, const Primitive& pr) const
    {
      return other_tree.do_intersect( internal::Primitive_helper<BaseTraits>::get_datum(pr,m_traits));
    }
    
    template<typename AABBTraits>
    bool operator()(const CGAL::AABB_tree<AABBTraits>& other_tree, const Bounding_box& bbox) const
    {
      Point_3 min(bbox.xmin(), bbox.ymin(), bbox.zmin()),
      max(bbox.xmax(), bbox.ymax(), bbox.zmax());
      
      min = m_traits.transformation().transform(min);
      max = m_traits.transformation().transform(max);
      Bounding_box transfo_box(to_double(min.x()), to_double(min.y()), to_double(min.z()),
                               to_double(max.x()), to_double(max.y()), to_double(max.z()));
      return other_tree.do_intersect(transfo_box);
    }
  };
  
  typedef typename BaseTraits::Intersection Intersection;
  
  //Distance Queries
  typedef typename BaseTraits::Compare_distance              Compare_distance;
  typedef typename BaseTraits::Closest_point                 Closest_point   ;
  typedef typename BaseTraits::Squared_distance              Squared_distance;
  typedef typename BaseTraits::Equal_3                       Equal_3         ;
  
  //Operations   
  Split_primitives split_primitives_object() const {
    return BaseTraits::split_primitives_object();
  }
  
  Compute_bbox compute_bbox_object() const{
    return BaseTraits::compute_bbox_object();
  }
  Do_intersect do_intersect_object() const{
    return Do_intersect(*this);
  }
  
  Intersection intersection_object() const{
    return BaseTraits::intersection_object();
  }
  
  Compare_distance compare_distance_object() const{
    return BaseTraits::compare_distance_object();
  }
  Closest_point closest_point_object() const{
    return BaseTraits::closest_point_object();
  }
  Squared_distance squared_distance_object() const{
    return BaseTraits::squared_distance_object();
  }
  Equal_3 equal_3_object() const{
    return BaseTraits::equal_3_object;
  }
  
  //Specific
  void set_transformation(const Aff_transformation_3<Kernel>& trans) const 
  {
    m_transfo = trans;
  }
  
  const Aff_transformation_3<Kernel>& transformation() const { return m_transfo; }
  
private:
  mutable Aff_transformation_3<Kernel> m_transfo;

};
  

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_TRANSFORMED_TRAITS_BASE_H
