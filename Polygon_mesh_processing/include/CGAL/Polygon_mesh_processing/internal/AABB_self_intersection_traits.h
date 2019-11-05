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

  //just for the declaration of the facet_intersector, it only needs to compile,
  // bc we use the operator() with faces, so the boxes are never used.
  template< typename Info>
  struct Dummy_box {
    Dummy_box(){}

    Info info() const
    {
      return Info();
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
      Dummy_box<typename boost::graph_traits<TM>::face_descriptor>, Output_iterator, VPM >
      (m, m_out_it, vpmap, Kernel());
  }
  ~AABB_self_intersection_traits() { delete facets_intersector;}

  bool go_further() const { return true; }

  void intersection(const Query& query, const Primitive& primitive)
  {
      if ( query.id() >= primitive.id()) return;
      facets_intersector->operator ()(query.id(), primitive.id());
  }

  bool do_intersect(const Query& query, const Node& node) const
  {

    return m_traits.do_intersect_object()(CGAL::internal::Primitive_helper<AABBTraits>::
                                          get_datum(query, m_traits), node.bbox());
  }

private:
  Output_iterator m_out_it;
  const AABBTraits& m_traits;
  CGAL::internal::Intersect_facets<TM, Kernel,
  Dummy_box<typename boost::graph_traits<TM>::face_descriptor>,
  Output_iterator, VPM > *facets_intersector;
};

//todo: move it in self_intersections.h, but careful to the Facet_intersectors that is in there and needed by this file.
// It makes it impossible to put that function in s_i.h without creating a dep-cycle right now. It probably needs to be moved in its own header.
template <class Concurrency_tag,
          typename Query, typename TM, typename VPM, typename Kernel,
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

#if !defined(CGAL_LINKED_WITH_TBB)
  CGAL_static_assertion_msg (!(boost::is_convertible<Concurrency_tag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if (boost::is_convertible<Concurrency_tag,Parallel_tag>::value)
  {
    std::vector<typename boost::graph_traits<TM>::face_descriptor> fs;
    fs.reserve(num_faces(m));
    for(const auto& f : faces(m))
    {
      fs.push_back(f);
    }
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0,num_faces(m)), [fs, m, traversal_traits, tree](const blocked_range<size_t>& r){
      for( size_t i=r.begin(); i!=r.end(); ++i )
      {
        Query q(fs[i], m);
        tree->traversal(q, traversal_traits);
      }
    });
  }
  else
#endif
  {
    for(const auto& f : faces(m))
    {
      Query q(f, m);
      tree->traversal(q, traversal_traits);
    }
  }
  return out;
}
}//end CGAL
#endif // CGAL_AABB_SELF_INTERSECTION_TRAITS_H
