// Copyright (c) 2016 GeometryFactory (France).
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
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_PREDICATES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_PREDICATES_H

#include <boost/graph/graph_traits.hpp>
#include <CGAL/use.h>

namespace CGAL {
namespace Corefinement {


template<class TriangleMesh, class VertexPointMap, class Nodes_vector>
struct Less_along_a_halfedge{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor hedge;
  const TriangleMesh& tm;
  const VertexPointMap& vpm;
  const Nodes_vector& nodes;

  Less_along_a_halfedge(halfedge_descriptor hedge_,
                        const TriangleMesh& tm_,
                        const VertexPointMap& vpm_,
                        const Nodes_vector& nodes_)
    : hedge(hedge_)
    , tm(tm_)
    , vpm(vpm_)
    , nodes(nodes_)
  {}

  bool operator()(std::size_t i, std::size_t j) const {
    //returns true, iff q lies strictly between p and r.
    typename Nodes_vector::Protector p;
    try{
      CGAL_USE(p);

      return CGAL::collinear_are_strictly_ordered_along_line(
        nodes.to_interval(get(vpm, target(hedge,tm))),
        nodes.interval_node(j),
        nodes.interval_node(i));
    }
    catch(CGAL::Uncertain_conversion_exception&){
      return CGAL::collinear_are_strictly_ordered_along_line(
        nodes.to_exact(get(vpm, target(hedge,tm))),
        nodes.exact_node(j),
        nodes.exact_node(i));
    }
  }
};

} } // end of namespace CGAL::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_PREDICATES_H
