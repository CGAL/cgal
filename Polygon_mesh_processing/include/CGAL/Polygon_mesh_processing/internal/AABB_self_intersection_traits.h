// Copyright (c) 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author        : Maxime Gimeno
//

#ifndef CGAL_AABB_SELF_INTERSECTION_TRAITS_H
#define CGAL_AABB_SELF_INTERSECTION_TRAITS_H

#include <CGAL/license/AABB_tree.h>


#include <CGAL/internal/AABB_tree/AABB_node.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h> // todo: only for the Facet Intersector. To be replaced if it moves.

namespace CGAL{


template<typename AABBTraits,class TM, class VPM, typename Kernel,
         typename Query, typename Output_iterator>
class AABB_self_intersection_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef CGAL::AABB_node<AABBTraits> Node;

  struct Dummy_box {

    Primitive_id id;

    Dummy_box(const Primitive_id& id)
      :id(id){}

    Primitive_id info() const
    {
      return id;
    }
  };

public:
  AABB_self_intersection_traits(const TM& m,
              VPM vpmap,
              Output_iterator out_it,
              const AABBTraits& traits)
    : m_out_it(out_it), m_traits(traits) {
  facets_intersector =
      new CGAL::internal::Intersect_facets<TM, typename AABBTraits::Geom_traits,
      Dummy_box, Output_iterator, VPM >
      (m, m_out_it, vpmap, Kernel());
  }
  ~AABB_self_intersection_traits() { delete facets_intersector;}

  bool go_further() const { return true; }

  void intersection(const Query& query, const Primitive& primitive)
  {
      if ( query.id() >= primitive.id()) return;
      Dummy_box a(query.id()), b(primitive.id());
      facets_intersector->operator ()(&a, &b);
  }

  bool do_intersect(const Query& query, const Node& node) const
  {

    return m_traits.do_intersect_object()(CGAL::internal::Primitive_helper<AABBTraits>::
                                          get_datum(query, m_traits), node.bbox());
  }

private:
  Output_iterator m_out_it;
  const AABBTraits& m_traits;
  CGAL::internal::Intersect_facets<TM, Kernel, Dummy_box, Output_iterator, VPM > *facets_intersector;
};

//todo: move it in self_intersections.h, but careful to the Facet_intersectors that is in there and needed by this file.
// It makes it impossible to put that function in s_i.h without creating a dep-cycle right now. It probably needs to be moved in its own header.
template<typename Query, typename TM, typename VPM, typename Kernel,
         typename Tr, typename OutputIterator> //query = face
OutputIterator self_intersections_with_tree(const TM& m,
                                  VPM vpmap,
                                  Kernel,
                                  CGAL::AABB_tree<Tr>* tree, //should probably be a template param to avoid include AABB_tree.h
                                  OutputIterator out) //const
{
  using namespace CGAL;

  AABB_self_intersection_traits<Tr,
      TM, VPM, Kernel,
      Query, OutputIterator> traversal_traits(m, vpmap, out,tree->traits());
  for(const auto& f : faces(m))
  {
    Query q(f, m);
    tree->traversal(q, traversal_traits);
  }
  return out;
}
}//end CGAL
#endif // CGAL_AABB_SELF_INTERSECTION_TRAITS_H
