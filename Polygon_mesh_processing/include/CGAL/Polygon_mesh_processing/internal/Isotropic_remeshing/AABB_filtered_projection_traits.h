// Copyright (c) 2008-2009  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// Copyright (c) 2010  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Camille Wormser, Pierre Alliez, Stephane Tayeb, Laurent Rineau
//
// File adapted from <CGAL/internal/AABB_tree/AABB_traversal_traits.h>
//

#ifndef CGAL_AABB_FILTERED_PROJECTION_TRAITS_H
#define CGAL_AABB_FILTERED_PROJECTION_TRAITS_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>


#include <CGAL/property_map.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/internal/AABB_tree/AABB_node.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

/**
 * @class Projection_traits
 */
template <typename AABBTraits,
          typename IndexPropertyMap,
          bool keep = false>
class Filtered_projection_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point_3;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;
  typedef typename ::CGAL::AABB_tree<AABBTraits>::size_type size_type;

  typedef typename boost::property_traits<IndexPropertyMap>::value_type Index_type;

  typedef std::set<typename boost::remove_const<Index_type>::type> Set_of_indices;

public:
  template <typename IndexToIgnoreIterator>
  Filtered_projection_traits(IndexToIgnoreIterator begin,
                             IndexToIgnoreIterator end,
                             const AABBTraits& aabb_traits,
                             IndexPropertyMap index_map = IndexPropertyMap())
    : m_closest_point_initialized(false),
      set_of_indices(begin, end),
      aabb_traits(aabb_traits),
      index_map(index_map)
  {
  }

  Filtered_projection_traits(Index_type index,
                             const AABBTraits& aabb_traits,
                             IndexPropertyMap index_map = IndexPropertyMap())
    : m_closest_point_initialized(false),
      set_of_indices(),
      aabb_traits(aabb_traits),
      index_map(index_map)
  {
    set_of_indices.insert(index);
  }

  bool go_further() const { return true; }

  void intersection(const Point_3& query, const Primitive& primitive)
  {
    const Index_type& id = get(index_map, primitive.id());

    if(keep != (set_of_indices.count(id) > 0)) return;
    if(!m_closest_point_initialized) {
      typedef CGAL::internal::Primitive_helper<AABBTraits> Helper;
      m_closest_point = Helper::get_reference_point(primitive, aabb_traits);
      m_closest_primitive = primitive.id();
      m_closest_point_initialized = true;
    }
    Point_3 new_closest_point = aabb_traits.closest_point_object()
      (query, primitive, m_closest_point);
    if(new_closest_point != m_closest_point)
    {
      m_closest_primitive = primitive.id();
      m_closest_point = new_closest_point; // this effectively shrinks the sphere
    }
  }

  bool do_intersect(const Point_3& query, const Node& node) const
  {
    if(!m_closest_point_initialized) return true;
    return AABBTraits().compare_distance_object()
      (query, node.bbox(), m_closest_point) == CGAL::SMALLER;
  }

  bool found() { return m_closest_point_initialized; };

  Point_3 closest_point() const { return m_closest_point; }
  Point_and_primitive_id closest_point_and_primitive() const
  {
    return Point_and_primitive_id(m_closest_point, m_closest_primitive);
  }

private:
  Point_3 m_closest_point;
  typename Primitive::Id m_closest_primitive;
  bool m_closest_point_initialized;
  Set_of_indices set_of_indices;
  const AABBTraits& aabb_traits;
  IndexPropertyMap index_map;
}; // end Filtered_projection_traits

} // end namespace internal
} // end namespace Polygon_mesh_processing
} // end namespace CGAL

#endif // CGAL_MESH_3_AABB_FILTERED_PROJECTION_TRAITS_H
