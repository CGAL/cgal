// Copyright (c) 2019 GeometryFactory (France).
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
// Author(s) : Sebastien Loriot, Martin Skrodzki
//

#ifndef CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE
#define CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE

#include <iostream>

namespace CGAL {

  /**
   * @class Hausdorff_primitive_traits
   */
  template<typename AABBTraits, typename Query>
  class Hausdorff_primitive_traits
  {
    typedef typename AABBTraits::Primitive Primitive;
    typedef ::CGAL::AABB_node<AABBTraits> Node;

  public:
    Hausdorff_primitive_traits(const AABBTraits& traits)
      : m_traits(traits) {}

    // Explore the whole tree, i.e. always enter children if the methods
    // do_intersect() below determines that it is worthwhile.
    bool go_further() const { return true; }

    // Compute the explicit Hausdorff distance to the given primitive
    void intersection(const Query& query, const Primitive& primitive)
    {
      // Have reached a single triangle
      std::cout << "Reached Triangle " << primitive.id() << '\n';

      /* TODO implement handling of a single triangle
      /  - Call Culling on B (First maybe don't cull, but only consider closest
      /    triangle), obtain local bounds for the triangle
      /  - Update global Hausdorff bounds according to the obtained local bounds
      /  - return the current best known global bounds
      */
    }

    bool do_intersect(const Query& query, const Node& node) const
    {
      // Have reached a node, determine whether or not to enter it

      /* TODO implement processing of an AABB node
      /  - Determine distance of the node's bounding box (node.bbox()) to the
      /    closest point in B (First maybe any point)
      /  - If the distance is larger than the global lower bound, enter the
      /    node, i.e. return true.
      */
      return true;
      //return m_traits.do_intersect_object()(query, node.bbox());
    }

  private:
    const AABBTraits& m_traits;
  };
}

#endif //CGAL_PMP_INTERNAL_AABB_TRAVERSAL_TRAITS_WITH_HAUSDORFF_DISTANCE
