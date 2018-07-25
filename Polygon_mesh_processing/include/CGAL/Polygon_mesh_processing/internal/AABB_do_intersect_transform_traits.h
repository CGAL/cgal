
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

#ifndef CGAL_AABB_DO_INTERSECT_TRANSFORM_TRAITS_H
#define CGAL_AABB_DO_INTERSECT_TRANSFORM_TRAITS_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <CGAL/internal/AABB_tree/Is_ray_intersection_geomtraits.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>
#include <CGAL/Filtered_predicate.h>

#include <CGAL/Aff_transformation_3.h>
#include <boost/mpl/if.hpp>

/// \file AABB_do_intersect_transform_traits.h

namespace CGAL {

template<typename Kernel, typename AABBPrimitive, 
         typename SUPPORTS_ROTATION = CGAL::Tag_true>
class AABB_do_intersect_transform_traits:
    public AABB_traits<Kernel, AABBPrimitive>
{
  mutable Aff_transformation_3<Kernel> m_transfo;
  mutable bool has_rotation;
  typedef AABB_traits<Kernel, AABBPrimitive> BaseTraits;
  typedef AABB_do_intersect_transform_traits<Kernel, AABBPrimitive, SUPPORTS_ROTATION> Self;
  
  void set_transformation(const Aff_transformation_3<Kernel>& trans, CGAL::Tag_true) const
  {
    has_rotation = (trans.m(0,1) != 0
        || trans.m(0,2) != 0
        || trans.m(1,0) != 0
        || trans.m(1,2) != 0
        || trans.m(2,0) != 0
        || trans.m(2,1) !=0);
  }
  void set_transformation(const Aff_transformation_3<Kernel>& trans, CGAL::Tag_false) const
  {}
  Bbox_3
  compute_transformed_bbox_impl(const Bbox_3& bbox, const Aff_transformation_3<Kernel>& transfo, Tag_true)const
  {
    if(has_rotation)
      return compute_transformed_bbox_impl(bbox, m_transfo, Tag_false());
    typedef Simple_cartesian<Interval_nt_advanced> AK;
    typedef Cartesian_converter<Kernel, AK>    C2F;
    C2F c2f;

    AK::Aff_transformation_3 af = c2f(transfo);

    AK::FT xtrm[6] = { c2f(bbox.min(0)), c2f(bbox.max(0)),
                       c2f(bbox.min(1)), c2f(bbox.max(1)),
                       c2f(bbox.min(2)), c2f(bbox.max(2)) };

    typename AK::Point_3 ps[8];
    ps[0] = af( AK::Point_3(xtrm[0], xtrm[2], xtrm[4]) );
    ps[1] = af( AK::Point_3(xtrm[0], xtrm[2], xtrm[5]) );
    ps[2] = af( AK::Point_3(xtrm[0], xtrm[3], xtrm[4]) );
    ps[3] = af( AK::Point_3(xtrm[0], xtrm[3], xtrm[5]) );

    ps[4] = af( AK::Point_3(xtrm[1], xtrm[2], xtrm[4]) );
    ps[5] = af( AK::Point_3(xtrm[1], xtrm[2], xtrm[5]) );
    ps[6] = af( AK::Point_3(xtrm[1], xtrm[3], xtrm[4]) );
    ps[7] = af( AK::Point_3(xtrm[1], xtrm[3], xtrm[5]) );

    return bbox_3(ps, ps+8);
  }

  Bbox_3
  compute_transformed_bbox_impl(const Bbox_3& bbox, const Aff_transformation_3<Kernel>& transfo, Tag_false)const
  {
    typedef Simple_cartesian<Interval_nt_advanced > AK;
    typedef Cartesian_converter<Kernel, AK>    C2F;
    C2F c2f;

    AK::Aff_transformation_3 af = c2f(transfo);

    AK::FT xtrm[6] = { c2f(bbox.min(0)), c2f(bbox.max(0)),
                       c2f(bbox.min(1)), c2f(bbox.max(1)),
                       c2f(bbox.min(2)), c2f(bbox.max(2)) };

    typename AK::Point_3 ps[2];
    ps[0] = af( AK::Point_3(xtrm[0], xtrm[2], xtrm[4]) );
    ps[1] = af( AK::Point_3(xtrm[1], xtrm[3], xtrm[5]) );

    return bbox_3(ps, ps+2);
  }
public:

  //Constructor
  AABB_do_intersect_transform_traits(const Aff_transformation_3<Kernel>& transf = Aff_transformation_3<Kernel>(IDENTITY))
  {
    has_rotation = false;
    set_transformation(transf, SUPPORTS_ROTATION());
  }

  // AABBTraits concept types
  typedef typename BaseTraits::Point_3 Point_3;
  typedef typename BaseTraits::Primitive Primitive;
  typedef typename BaseTraits::Bounding_box Bounding_box;

  // helper functions

  
  Bbox_3
  compute_transformed_bbox(const Bbox_3& bbox) const
  {
    return compute_transformed_bbox_impl(bbox, m_transfo, SUPPORTS_ROTATION());
  }

  // Do_intersect predicate
  class Do_intersect
    : BaseTraits::Do_intersect
  {
    typedef AABB_do_intersect_transform_traits<Kernel, AABBPrimitive, SUPPORTS_ROTATION> AABBTraits;
    const AABBTraits& m_traits;
    typedef typename BaseTraits::Do_intersect Base;

    Bounding_box
    compute_transformed_bbox(const Bounding_box& bbox) const
    {
      return m_traits.compute_transformed_bbox(bbox);
    }

  public:
    Do_intersect(const AABBTraits& traits)
    : Base(static_cast<const BaseTraits&>(traits)),
      m_traits(traits)
    {}

    template<typename Query>
    bool operator()(const Query& q, const Bounding_box& bbox) const
    {
      return
        static_cast<const Base*>(this)->operator()(
          q, compute_transformed_bbox(bbox));
    }

    template<typename Query>
    bool operator()(const Query& q, const Primitive& pr) const
    {
      // transformation is done within Primitive_helper
      return do_intersect(q, internal::Primitive_helper<Self>::get_datum(pr,m_traits));
    }

    // intersection with AABB-tree
    template<typename AABBTraits>
    bool operator()(const CGAL::AABB_tree<AABBTraits>& other_tree, const Primitive& pr) const
    {
      // transformation is done within Primitive_helper
      return other_tree.do_intersect( internal::Primitive_helper<Self>::get_datum(pr,m_traits));
    }

    template<typename AABBTraits>
    bool operator()(const CGAL::AABB_tree<AABBTraits>& other_tree, const Bounding_box& bbox) const
    {
      return other_tree.do_intersect(compute_transformed_bbox(bbox));
    }
  };

  Do_intersect do_intersect_object() const{
    return Do_intersect(*this);
  }

  //Specific
  void set_transformation(const Aff_transformation_3<Kernel>& trans) const
  {
    m_transfo = trans;
    set_transformation(trans, SUPPORTS_ROTATION());
  }

  const Aff_transformation_3<Kernel>& transformation() const { return m_transfo; }
};

namespace internal {

template<typename K, typename P, typename T>
struct Primitive_helper<AABB_do_intersect_transform_traits<K,P,T> ,true>{

typedef AABB_do_intersect_transform_traits<K,P,T> Traits;


static typename Traits::Primitive::Datum get_datum(const typename Traits::Primitive& p,
                            const Traits & traits)
{
  return p.datum(traits.shared_data()).transform(traits.transformation());
}

static typename Traits::Point_3 get_reference_point(const typename Traits::Primitive& p,const Traits& traits) {
  return p.reference_point(traits.shared_data()).transform(traits.transformation());
}

};

template<typename K, typename P, typename T>
typename CGAL::AABB_tree<AABB_do_intersect_transform_traits<K,P,T> >::Bounding_box
get_tree_bbox(const CGAL::AABB_tree<AABB_do_intersect_transform_traits<K,P,T> >& tree)
{
  return tree.traits().compute_transformed_bbox(tree.bbox());
}

template<typename K, typename P, typename T>
typename CGAL::AABB_tree<AABB_do_intersect_transform_traits<K,P,T> >::Bounding_box
get_node_bbox(const CGAL::AABB_node<AABB_do_intersect_transform_traits<K,P,T> >& node,
              const AABB_do_intersect_transform_traits<K,P,T>& traits)
{
  return traits.compute_transformed_bbox(node.bbox());
}

} // end internal

}//end CGAL

#endif //CGAL_AABB_AABB_do_intersect_transform_traits_H
