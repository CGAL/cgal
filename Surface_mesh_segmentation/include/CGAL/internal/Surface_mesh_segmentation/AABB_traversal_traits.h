// Copyright (c) 2014  GeometryFactory Sarl (France).
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
// Author(s)     : Ilker O. Yaz


#ifndef CGAL_INTERNAL_SURFACE_MESH_SEGMENTATION_AABB_TRAVERSAL_TRAITS_H
#define CGAL_INTERNAL_SURFACE_MESH_SEGMENTATION_AABB_TRAVERSAL_TRAITS_H

#include <CGAL/license/Surface_mesh_segmentation.h>


namespace CGAL
{

/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @class Special case for ray/segment-triangle
 * the only difference with the offical one (Listing_intersection_traits) is that
 * is the do_intersect which is made prior to the intersection call.
 */
template<typename AABBTraits, typename Query, typename Output_iterator>
class Listing_intersection_traits_ray_or_segment_triangle
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;
  typedef typename ::CGAL::AABB_tree<AABBTraits>::size_type size_type;

public:
  Listing_intersection_traits_ray_or_segment_triangle(Output_iterator out_it,
      const AABBTraits& traits)
    : m_out_it(out_it), m_traits(traits) {}

  bool go_further() const {
    return true;
  }

  void intersection(const Query& query, const Primitive& primitive) {
    //SL: using Kernel_traits is not bad in this context cause we expect a Ray/Segment from a CGAL Kernel here
    typedef typename Kernel_traits<Query>::Kernel GeomTraits;
    typedef typename AABBTraits:: template Intersection_and_primitive_id<Query>::Type Intersection_and_primitive_id;

    if ( GeomTraits().do_intersect_3_object()(query,
         primitive.datum(m_traits.shared_data())) ) {
      boost::optional<Intersection_and_primitive_id> intersection
        = m_traits.intersection_object()(query, primitive);
      if(intersection) {
        *m_out_it++ = *intersection;
      }
    }
  }

  bool do_intersect(const Query& query, const Node& node) const {
    return m_traits.do_intersect_object()(query, node.bbox());
  }

private:
  Output_iterator m_out_it;
  const AABBTraits& m_traits;
};

/// @endcond

} //namespace CGAL
#endif
