// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Maxime Gimeno
//             Sebastien Loriot
//

#ifndef CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_TRANSFORMATION
#define CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_TRANSFORMATION

#include <CGAL/license/Polygon_mesh_processing/collision_detection.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <CGAL/internal/AABB_tree/Is_ray_intersection_geomtraits.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Aff_transformation_3.h>
#include <boost/mpl/if.hpp>

namespace CGAL {

template<typename AABBTraits, typename Kernel, class SUPPORTS_ROTATION>
class Traversal_traits_with_transformation_helper
{
  Bbox_3
  compute_transformed_bbox_impl(const CGAL::Aff_transformation_3<Kernel>& at,
                                const Bbox_3& bbox,
                                bool has_rotation,
                                /*SUPPORTS_ROTATION*/ Tag_true) const
  {
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_UPWARD);

    if(!has_rotation)
      return compute_transformed_bbox_impl(at, bbox, has_rotation, Tag_false());

    typedef Simple_cartesian<Interval_nt_advanced> AK;
    typedef Cartesian_converter<Kernel, AK>    C2F;
    C2F c2f;

    AK::Aff_transformation_3 a_at = c2f(at);

    AK::FT xtrm[6] = { c2f((bbox.min)(0)), c2f((bbox.max)(0)),
                       c2f((bbox.min)(1)), c2f((bbox.max)(1)),
                       c2f((bbox.min)(2)), c2f((bbox.max)(2)) };

    typename AK::Point_3 ps[8];
    ps[0] = a_at( AK::Point_3(xtrm[0], xtrm[2], xtrm[4]) );
    ps[1] = a_at( AK::Point_3(xtrm[0], xtrm[2], xtrm[5]) );
    ps[2] = a_at( AK::Point_3(xtrm[0], xtrm[3], xtrm[4]) );
    ps[3] = a_at( AK::Point_3(xtrm[0], xtrm[3], xtrm[5]) );

    ps[4] = a_at( AK::Point_3(xtrm[1], xtrm[2], xtrm[4]) );
    ps[5] = a_at( AK::Point_3(xtrm[1], xtrm[2], xtrm[5]) );
    ps[6] = a_at( AK::Point_3(xtrm[1], xtrm[3], xtrm[4]) );
    ps[7] = a_at( AK::Point_3(xtrm[1], xtrm[3], xtrm[5]) );

    return bbox_3(ps, ps+8);
  }

  Bbox_3
  compute_transformed_bbox_impl(const CGAL::Aff_transformation_3<Kernel>& at,
                                const Bbox_3& bbox,
                                bool,
                                /*SUPPORTS_ROTATION*/ Tag_false) const
  {
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_UPWARD);

    typedef Simple_cartesian<Interval_nt_advanced > AK;
    typedef Cartesian_converter<Kernel, AK>    C2F;
    C2F c2f;

    AK::Aff_transformation_3 a_at = c2f(at);

    AK::FT xtrm[6] = { c2f((bbox.min)(0)), c2f((bbox.max)(0)),
                       c2f((bbox.min)(1)), c2f((bbox.max)(1)),
                       c2f((bbox.min)(2)), c2f((bbox.max)(2)) };

    typename AK::Point_3 ps[2];
    ps[0] = a_at( AK::Point_3(xtrm[0], xtrm[2], xtrm[4]) );
    ps[1] = a_at( AK::Point_3(xtrm[1], xtrm[3], xtrm[5]) );

    return bbox_3(ps, ps+2);
  }

public:

  bool has_rotation(const CGAL::Aff_transformation_3<Kernel>& at) const
  {
    return  ( at.m(0,1) != 0 || at.m(0,2) != 0 || at.m(1,0) != 0
           || at.m(1,2) != 0 || at.m(2,0) != 0 || at.m(2,1) !=0);
  }

  Bbox_3
  compute_transformed_bbox(const CGAL::Aff_transformation_3<Kernel>& at,
                           const Bbox_3& bbox,
                           bool has_rotation) const
  {
    return compute_transformed_bbox_impl(at, bbox, has_rotation, SUPPORTS_ROTATION());
  }
};

// traversal traits for a tree vs a primitive
template<typename AABBTraits, typename Kernel,
         typename SUPPORTS_ROTATION = CGAL::Tag_true>
class Do_intersect_traversal_traits_with_transformation
  : public Traversal_traits_with_transformation_helper<AABBTraits, Kernel, SUPPORTS_ROTATION>
{
  typedef typename AABBTraits::Primitive Primitive;
  typedef ::CGAL::AABB_node<AABBTraits> Node;
  typedef Traversal_traits_with_transformation_helper<AABBTraits, Kernel, SUPPORTS_ROTATION> Base;

  void register_transformation(CGAL::Tag_true)
  {
    m_has_rotation = this->has_rotation(m_transfo);
  }

  void register_transformation(CGAL::Tag_false)
  {}

public:
  Do_intersect_traversal_traits_with_transformation():
    m_traits_ptr(nullptr)
  {}

  Do_intersect_traversal_traits_with_transformation(const AABBTraits& traits)
    : m_is_found(false), m_traits_ptr(&traits), m_has_rotation(false)
  {}

  bool go_further() const { return !m_is_found; }

  template <class Query>
  void intersection(const Query& query, const Primitive& primitive)
  {
    if( CGAL::do_intersect(query,
                     internal::Primitive_helper<AABBTraits>::get_datum(primitive, *m_traits_ptr).transform(m_transfo)) )
      m_is_found = true;
  }

  template <class Query>
  bool do_intersect(const Query& query, const Node& node) const
  {
    return m_traits_ptr->do_intersect_object()(query, compute_transformed_bbox(node.bbox()));
  }

  bool is_intersection_found() const { return m_is_found; }

  void reset()
  {
    m_is_found = false;
  }

  const Aff_transformation_3<Kernel>&
  transformation() const
  {
    return m_transfo;
  }

  void set_transformation(const Aff_transformation_3<Kernel>& transfo)
  {
    m_transfo = transfo;
    register_transformation(SUPPORTS_ROTATION());
  }

  Bbox_3
  compute_transformed_bbox(const Bbox_3& bbox) const
  {
    return Base::compute_transformed_bbox(m_transfo, bbox, m_has_rotation);
  }

  // helper for Point_inside_vertical_ray_cast
  class Transformed_tree_helper
  {
    typedef AABB_tree<AABBTraits> Tree;
    typedef CGAL::AABB_node<AABBTraits> Node;
    typedef Do_intersect_traversal_traits_with_transformation<AABBTraits, Kernel, SUPPORTS_ROTATION> Traversal_traits;

    Traversal_traits m_tt;

  public:

    Transformed_tree_helper(const Traversal_traits& tt)
      : m_tt(tt)
    {}

    Bbox_3 get_tree_bbox(const AABB_tree<AABBTraits>& tree) const
    {
      return m_tt.compute_transformed_bbox(tree.bbox());
    }

    typename AABBTraits::Primitive::Datum
    get_primitive_datum(const typename AABBTraits::Primitive& primitive, const AABBTraits& traits) const
    {
      return internal::Primitive_helper<AABBTraits>::get_datum(primitive, traits).transform(m_tt.transformation());
    }

    Bbox_3 get_node_bbox(const Node& node) const
    {
      return m_tt.compute_transformed_bbox(node.bbox());
    }
  };

  Transformed_tree_helper get_helper() const
  {
    return Transformed_tree_helper(*this);
  }

private:
  bool m_is_found;
  const AABBTraits* m_traits_ptr;
  Aff_transformation_3<Kernel> m_transfo;
  bool m_has_rotation;
};


// traversal traits for a tree
template<typename AABBTraits, class Kernel, typename SUPPORTS_ROTATION = CGAL::Tag_true >
class Do_intersect_traversal_traits_for_two_trees
  : public Traversal_traits_with_transformation_helper<AABBTraits, Kernel, SUPPORTS_ROTATION>
{
  typedef typename AABBTraits::Primitive Primitive;
  typedef ::CGAL::AABB_node<AABBTraits> Node;
  typedef Traversal_traits_with_transformation_helper<AABBTraits, Kernel, SUPPORTS_ROTATION> Base;
  typedef Do_intersect_traversal_traits_with_transformation<AABBTraits, Kernel, SUPPORTS_ROTATION> Query_traversal_traits;

  void register_transformation(CGAL::Tag_true)
  {
    m_has_rotation = this->has_rotation(m_transfo);
  }

  void register_transformation(CGAL::Tag_false)
  {}

  Bbox_3
  compute_transformed_bbox(const Bbox_3& bbox) const
  {
    return Base::compute_transformed_bbox(m_transfo, bbox, m_has_rotation);
  }

public:
  Do_intersect_traversal_traits_for_two_trees(const AABBTraits& traits,
                                    const Aff_transformation_3<Kernel>& transfo,
                                    const Query_traversal_traits& query_traversal_traits)
    : m_is_found(false)
    , m_traits(traits)
    , m_transfo(transfo)
    , m_has_rotation(false)
    , m_query_traversal_traits(query_traversal_traits)

  {
    register_transformation(SUPPORTS_ROTATION());
  }

  bool go_further() const { return !m_is_found; }

  void intersection(const AABB_tree<AABBTraits>& query, const Primitive& primitive)
  {
    query.traversal( internal::Primitive_helper<AABBTraits>::get_datum(primitive,m_traits).transform(m_transfo), m_query_traversal_traits );
    m_is_found = m_query_traversal_traits.is_intersection_found();
    m_query_traversal_traits.reset();
  }

  bool do_intersect(const AABB_tree<AABBTraits>& query, const Node& node)
  {
    query.traversal( compute_transformed_bbox(node.bbox()), m_query_traversal_traits );
    bool res = m_query_traversal_traits.is_intersection_found();
    m_query_traversal_traits.reset();
    return res;
  }

  bool is_intersection_found() const { return m_is_found; }

private:
  bool m_is_found;
  const AABBTraits& m_traits;
  const Aff_transformation_3<Kernel>& m_transfo;
  bool m_has_rotation;
  Do_intersect_traversal_traits_with_transformation<AABBTraits, Kernel, SUPPORTS_ROTATION> m_query_traversal_traits;
};

}//end CGAL

#endif //CGAL_AABB_AABB_do_intersect_transform_traits_H
