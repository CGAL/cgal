
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

/// \file AABB_transformed_traits.h

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
           typename Kernel,
           class Has_filtered_predicates = typename Kernel::Has_filtered_predicates_tag /* Tag_false*/>
class AABB_transformed_traits:
  public BaseTraits
{
public:

  //Constructor
  AABB_transformed_traits(const Aff_transformation_3<Kernel>& transf = Aff_transformation_3<Kernel>(IDENTITY))
    :m_transfo(transf)
  {}
  // AABBTraits concept types
  typedef typename BaseTraits::Point_3 Point_3;
  typedef typename BaseTraits::Primitive Primitive;
  typedef typename BaseTraits::Bounding_box Bounding_box;
  //Intersections
  class Do_intersect {
    const AABB_transformed_traits<BaseTraits, Kernel>& m_traits;
  public:
    Do_intersect(const AABB_transformed_traits<BaseTraits, Kernel>& traits)
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
  
  Do_intersect do_intersect_object() const{
    return Do_intersect(*this);
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
  
template<typename BaseTraits, 
         typename Kernel>
class AABB_transformed_traits<BaseTraits, Kernel, Tag_true>:
public BaseTraits
{
  typedef typename Kernel::Exact_kernel EK;
  typedef typename Kernel::Approximate_kernel AK;
  typedef typename Kernel::C2E C2E;
  typedef typename Kernel::C2F C2A;
public:

  //Constructor
  AABB_transformed_traits(const Aff_transformation_3<Kernel>& transf = Aff_transformation_3<Kernel>(IDENTITY))
    :m_transfo(transf)
  {}
  // AABBTraits concept types
  typedef typename BaseTraits::Point_3 Point_3;
  typedef typename BaseTraits::Primitive Primitive;
  typedef typename BaseTraits::Bounding_box Bounding_box;
  //Intersections
private:
  class Unfiltered_do_intersect {
    const AABB_transformed_traits<BaseTraits, Kernel>& m_traits;
  public:
    Unfiltered_do_intersect(const AABB_transformed_traits<BaseTraits, Kernel>& traits)
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
public:
  typedef Unfiltered_do_intersect Do_intersect;
  typedef typename BaseTraits::Intersection Intersection;
  
  Do_intersect do_intersect_object() const{
    return Do_intersect(*this);
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

#endif // CGAL_AABB_TRANSFORMED_TRAITS_H
