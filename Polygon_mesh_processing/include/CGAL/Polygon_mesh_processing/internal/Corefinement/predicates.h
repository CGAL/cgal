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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_PREDICATES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_PREDICATES_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <boost/graph/graph_traits.hpp>
#include <CGAL/use.h>

namespace CGAL {
namespace Corefinement {


template<class TriangleMesh, class VertexPointMap, class Node_vector>
struct Less_along_a_halfedge{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor hedge;
  const TriangleMesh& tm;
  const VertexPointMap& vpm;
  const Node_vector& nodes;

  Less_along_a_halfedge(halfedge_descriptor hedge_,
                        const TriangleMesh& tm_,
                        const VertexPointMap& vpm_,
                        const Node_vector& nodes_)
    : hedge(hedge_)
    , tm(tm_)
    , vpm(vpm_)
    , nodes(nodes_)
  {}

  bool operator()(std::size_t i, std::size_t j) const {
    return CGAL::collinear_are_strictly_ordered_along_line(
      nodes.to_exact(get(vpm, target(hedge,tm))),
      nodes.exact_node(j),
      nodes.exact_node(i));
  }
};

//Considering the plane with normal vector [o_prime,o] and containing o.
//We define the counterclockwise order around o when looking from
//the side of the plane containing o_prime.
//We consider the portion of the plane defined by rotating a ray starting at o
//from the planar projection of p1 to the planar projection of p2 in
//counterclockwise order.
//The predicates indicates whether the planar projection of point q lies in this
//portion of the plane.
//Preconditions:
//  o_prime,o,p1 are not collinear
//  o_prime,o,p2 are not collinear
//  o_prime,o,q are not collinear
//  o_prime,o,p1,q are not coplanar or coplanar_orientation(o,o_prime,p1,q)==NEGATIVE
//  o_prime,o,p2,q are not coplanar or coplanar_orientation(o,o_prime,p2,q)==NEGATIVE
template <class Kernel>
bool  sorted_around_edge(
  const typename Kernel::Point_3& o_prime, const typename Kernel::Point_3& o,
  const typename Kernel::Point_3& p1, const typename Kernel::Point_3& p2,
  const typename Kernel::Point_3& q)
{
  //guarantee to have non-flat triangles
  CGAL_precondition( !collinear(o_prime,o,p1) );
  CGAL_precondition( !collinear(o_prime,o,p2) );
  CGAL_precondition( !collinear(o_prime,o,q)  );

  //no two triangles are coplanar and on the same side of their common edge
  CGAL_precondition( !coplanar(o_prime,o,p1,q)
                     || coplanar_orientation(o,o_prime,p1,q)==NEGATIVE );
  CGAL_precondition( !coplanar(o_prime,o,p2,q)
                     || coplanar_orientation(o,o_prime,p2,q)==NEGATIVE );

  typename Kernel::Orientation_3 orientation;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename cpp11::result_of<
    typename Kernel::Orientation_3(Point_3, Point_3, Point_3, Point_3)>::type
      Orientation;

  Orientation s0 = orientation(o_prime, o, p1, p2);

  if ( s0==COPLANAR ) {
    //o, o_prime, p1 and p2 are coplanar
    Orientation orient=orientation(o_prime,o,p1,q);
    CGAL_precondition(orient!=COPLANAR);
    return orient==POSITIVE;
  }

  //o, o_prime, p1 and p2 are not coplanar
  Orientation s1 = orientation(o_prime, o, p1, q);
  Orientation s2 = orientation(o_prime, o, q , p2);

  if (s0 == POSITIVE) // the angle p1,o,p2 is smaller that Pi.
    return ( s1 == POSITIVE )
           && ( s2 ==POSITIVE ); //true if the angles p1,o,q and q,o,p2 are smaller than Pi
  else
    return ( s1 != NEGATIVE )
           || ( s2 != NEGATIVE ); //true if the angle p1,o,q or the angle q,o,p2 is smaller than or equal to Pi
}

template <class Kernel>
bool  are_triangles_coplanar_same_side(
    const typename Kernel::Point_3& o_prime, const typename Kernel::Point_3& o,
    const typename Kernel::Point_3& P, const typename Kernel::Point_3& q)
{
  if ( CGAL::orientation(o_prime, o, P ,q) != COPLANAR )
    return false;
  Orientation cpl_orient = coplanar_orientation(o_prime, o, P, q);
  CGAL_assertion( cpl_orient != COLLINEAR );
  return cpl_orient == POSITIVE;
}


template <class Node_id, class Node_vector, class vertex_descriptor, class Vpm>
bool are_triangles_coplanar_same_side(Node_id o_prime_index,
                                      Node_id o_index,
                                      Node_id p_index,
                                      Node_id q_index,
                                      vertex_descriptor p,
                                      vertex_descriptor q,
                                      const Vpm& vpm_p,
                                      const Vpm& vpm_q,
                                      const Node_vector& nodes)
{
  const Node_id NID((std::numeric_limits<Node_id>::max)());
    return are_triangles_coplanar_same_side<typename Node_vector::Exact_kernel>(
      nodes.exact_node(o_prime_index),
      nodes.exact_node(o_index),
      p_index == NID ? nodes.to_exact(get(vpm_p,p)) : nodes.exact_node(p_index),
      q_index == NID ? nodes.to_exact(get(vpm_q,q)) : nodes.exact_node(q_index)
    );
}

template <class Node_id, class Node_vector, class vertex_descriptor, class Vpm>
bool sorted_around_edge( Node_id o_prime_index,
                         Node_id o_index,
                         Node_id p1_index,
                         Node_id p2_index,
                         Node_id q_index,
                         vertex_descriptor p1,
                         vertex_descriptor p2,
                         vertex_descriptor q,
                         const Vpm& vpm_p,
                         const Vpm& vpm_q,
                         const Node_vector& nodes)
{
  const Node_id NID((std::numeric_limits<Node_id>::max)());
  return sorted_around_edge<typename Node_vector::Exact_kernel>(
           nodes.exact_node(o_prime_index),
           nodes.exact_node(o_index),
           p1_index == NID ? nodes.to_exact(get(vpm_p,p1))
                           : nodes.exact_node(p1_index),
           p2_index == NID ? nodes.to_exact(get(vpm_p,p2))
                           : nodes.exact_node(p2_index),
           q_index  == NID ? nodes.to_exact(get(vpm_q,q))
                           : nodes.exact_node(q_index ) );
}


} } // end of namespace CGAL::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_PREDICATES_H
